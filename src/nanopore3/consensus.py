"""Deterministic consensus construction with optional external-tool backends.

The default ``portable`` backend deliberately uses only edlib.  ``mafft_spoa`` is
an opt-in compatibility backend for the legacy Nanopore2 workflow; importing and
using the package never requires either external executable unless that backend is
selected explicitly.
"""

from __future__ import annotations

import hashlib
import math
import re
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import edlib

from .errors import Nanopore3Error

_CIGAR_TOKEN = re.compile(r"(\d+)([=XID])")
SUPPORTED_CONSENSUS_BACKENDS = frozenset({"portable", "mafft_spoa"})


class ConsensusBackendUnavailable(Nanopore3Error, RuntimeError):
    """Raised when an explicitly selected optional backend cannot be executed."""


class ConsensusBackendError(Nanopore3Error, RuntimeError):
    """Raised when an optional consensus tool fails or emits unusable output."""


@dataclass(frozen=True, slots=True)
class ConsensusRead:
    """A read eligible to contribute to one consensus."""

    read_uid: str
    sequence: str
    qualities: tuple[int, ...] | None = None


@dataclass(frozen=True, slots=True)
class ConsensusResult:
    sequence: str
    status: str
    n_reads_available: int
    n_reads_used: int
    contributor_ids: tuple[str, ...]
    mean_depth: float
    min_depth: int
    ambiguous_bases: int
    backend: str = "portable"
    failure_reason: str | None = None
    # Positions where more than one allele beat the background error rate, the
    # rate itself, the weakest winning support, and base observations dropped for
    # quality. All zero for backends that do not compute them.
    mixed_positions: int = 0
    # One entry per mixed position: (reference index, reference base, the
    # competing alleles sorted, their read fractions in the same order). Which
    # alleles disagreed distinguishes a damage signature from a genuine two-clone
    # well, and the fractions answer the question screening actually asks - is the
    # designed sequence still present in this well, and at what share. An "N" says
    # neither.
    mixed_alleles: tuple[tuple[int, str, str, tuple[float, ...]], ...] = ()
    background_error_rate: float = 0.0
    deletion_error_rate: float = 0.0
    weakest_support: float = 0.0
    low_quality_bases: int = 0


def _stable_rank(read_uid: str, seed: int, group_id: str) -> str:
    value = f"{seed}\0{group_id}\0{read_uid}".encode()
    return hashlib.sha256(value).hexdigest()


def select_reads(
    reads: Iterable[ConsensusRead], *, max_reads: int, seed: int, group_id: str
) -> tuple[ConsensusRead, ...]:
    """Select contributors independently of FASTQ order or worker scheduling."""

    unique: dict[str, ConsensusRead] = {}
    for read in reads:
        if read.read_uid in unique:
            raise ValueError(f"duplicate read_uid in consensus group: {read.read_uid}")
        unique[read.read_uid] = read
    ranked = sorted(
        unique.values(), key=lambda r: (_stable_rank(r.read_uid, seed, group_id), r.read_uid)
    )
    return tuple(ranked[:max_reads])


def validate_consensus_backend(backend: str) -> str:
    """Validate and return a public consensus backend name."""

    if backend not in SUPPORTED_CONSENSUS_BACKENDS:
        supported = ", ".join(sorted(SUPPORTED_CONSENSUS_BACKENDS))
        raise ValueError(
            f"unsupported consensus backend {backend!r}; choose one of: {supported}"
        )
    return backend


def _resolve_executable(name: str, configured: str | Path | None) -> str | None:
    candidate = str(configured) if configured is not None else name
    located = shutil.which(candidate)
    return str(Path(located).resolve()) if located is not None else None


def require_mafft_spoa(
    *, mafft_path: str | Path | None = None, spoa_path: str | Path | None = None
) -> tuple[str, str]:
    """Resolve both optional tools or raise one actionable installation error."""

    mafft = _resolve_executable("mafft", mafft_path)
    spoa = _resolve_executable("spoa", spoa_path)
    missing = [name for name, path in (("mafft", mafft), ("spoa", spoa)) if path is None]
    if missing:
        raise ConsensusBackendUnavailable(
            "consensus backend 'mafft_spoa' requires both MAFFT and SPOA on PATH; "
            f"missing: {', '.join(missing)}. Install the missing executable(s), "
            "select consensus.backend: portable, or provide explicit executable paths."
        )
    return mafft, spoa


# Estimated per-group error rates outside this range are not believable and would
# make every position either trivially significant or trivially not.
_MIN_ERROR_RATE, _MAX_ERROR_RATE = 0.002, 0.25


def binomial_upper_tail(successes: int, trials: int, probability: float) -> float:
    """P(X >= successes) for X ~ Binomial(trials, probability).

    Summed from the observed count upwards, which converges in a few terms in the
    upper tail where it is used.  Written out rather than taken from scipy so the
    portable backend keeps its only dependency being edlib.
    """

    if successes <= 0:
        return 1.0
    if successes > trials:
        return 0.0
    if probability <= 0.0:
        return 0.0
    if probability >= 1.0:
        return 1.0
    log_term = (
        math.lgamma(trials + 1)
        - math.lgamma(successes + 1)
        - math.lgamma(trials - successes + 1)
        + successes * math.log(probability)
        + (trials - successes) * math.log1p(-probability)
    )
    term = math.exp(log_term)
    total = term
    index = successes
    ratio = probability / (1.0 - probability)
    while index < trials:
        term *= (trials - index) / (index + 1) * ratio
        index += 1
        total += term
        if term <= total * 1e-15:
            break
    return min(1.0, total)


def _majority_call(
    counts: Counter[str],
    *,
    depth: int,
    min_support: float,
) -> tuple[str, bool, float]:
    """Plain majority, for polishing a draft rather than judging a clone.

    Iterative polishing wants the reads' own consensus with no reference bias and
    no abstentions: an ``N`` mid-refinement removes information the next round
    needs.  Significance testing belongs where a clone is being reported, not
    where a scaffold is being improved.
    """

    if not counts or depth <= 0:
        return "N", False, 0.0
    winner, winning = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
    support = winning / depth
    if support < min_support:
        return "N", False, support
    return winner, False, support


def _call_position(
    counts: Counter[str],
    reference_base: str,
    *,
    depth: int,
    spanning: int,
    substitution_error: float,
    deletion_error: float,
    positions: int,
    significance: float,
    minor_fraction: float,
    deletion_support: float,
) -> tuple[str, bool, float, tuple[str, tuple[float, ...]]]:
    """Decide one position: the base, whether it is mixed, its support, and alleles.

    The fourth value is the competing alleles when the position is mixed - sorted
    and joined, with their read fractions in the same order - and empty otherwise.

    An allele has to clear two independent bars: a one-sided binomial test against
    the group's own background error rate, and a minimum share of the reads.

    Both are needed.  A flat support fraction alone accepts 61% against 32% as a
    confident base, when that split is two alleles.  The test alone is worse: it
    assumes errors are independent at a uniform rate, and nanopore errors are
    neither - they cluster in homopolymers and awkward contexts, so at a hard
    position 15% of reads can disagree *systematically*.  Run on its own at 300x
    depth the test called nearly every position mixed.  The fraction floor is what
    makes it usable.

    Only substitutions are judged this way.  Measured on this data the mean
    substitution disagreement is 0.23% while the mean deletion disagreement is
    1.1%, and the deletion rate is wildly position-specific: individual positions
    run 20-25% deletions with no mixture present, which is the platform's error
    mode and not something a single per-group rate can describe.  Testing
    deletions against a scalar rate made 1,778 of 1,789 contested positions
    deletions.  So a deletion is taken only when it is the outright majority, as
    before, and never reported as mixture; a substitution minority at 20% is a
    hundred times its background and cannot be error.
    """

    if not counts or depth <= 0:
        return "N", False, 0.0, ("", ())

    # A deletion is a majority decision at the same threshold the majority caller
    # uses, never a mixture call. Taking one at >50% instead silently shortened 76
    # clones in a full run, because nanopore deletion rates sit in that band.
    deletions = counts.get("-", 0)
    bases = Counter({a: c for a, c in counts.items() if a != "-"})
    # Against every read that spanned the position, not only those whose base
    # survived the quality filter.
    reads_here = spanning or depth
    if deletions and deletions / reads_here >= deletion_support:
        return "-", False, deletions / reads_here, ("", ())
    if not bases:
        return "N", False, deletions / depth, ("", ())

    winner, winning_count = sorted(bases.items(), key=lambda item: (-item[1], item[0]))[0]
    support = winning_count / depth

    # Which substitutions are more common than substitution error explains?
    significant = [
        allele
        for allele, count in bases.items()
        if count >= 2
        and count / depth >= minor_fraction
        and count > substitution_error * depth
        and binomial_upper_tail(count, depth, substitution_error) * positions
        < significance
    ]
    if len(significant) > 1:
        ordered = sorted(significant)
        return (
            "N",
            True,
            support,
            ("".join(ordered), tuple(bases[allele] / depth for allele in ordered)),
        )
    if winner == reference_base:
        return winner, False, support, ("", ())
    # A non-reference winner has to be significant in its own right; otherwise the
    # position is noise and saying nothing is more honest than picking a side.
    if winner in significant:
        return winner, False, support, ("", ())
    return "N", False, support, ("", ())


def _median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def _build_portable_consensus(
    reference: str,
    reads: Sequence[ConsensusRead],
    *,
    group_id: str,
    min_depth: int = 3,
    max_reads: int = 100,
    min_support: float = 0.60,
    seed: int = 0,
    min_base_quality: int = 10,
    significance: float = 0.05,
    minor_fraction: float = 0.20,
    caller: str = "statistical",
) -> ConsensusResult:
    """Build the edlib reference-guided consensus used by the portable backend."""

    ref = reference.upper()
    if not ref:
        raise ValueError("reference cannot be empty")
    if min_depth < 1 or max_reads < 1:
        raise ValueError("min_depth and max_reads must be positive")
    if not 0.5 <= min_support <= 1.0:
        raise ValueError("min_support must be between 0.5 and 1.0")
    if not 0.0 < significance < 1.0:
        raise ValueError("significance must be between 0 and 1")

    available = len(reads)
    selected = select_reads(reads, max_reads=max_reads, seed=seed, group_id=group_id)
    if len(selected) < min_depth:
        return ConsensusResult(
            sequence="",
            status="low_depth",
            n_reads_available=available,
            n_reads_used=len(selected),
            contributor_ids=tuple(r.read_uid for r in selected),
            mean_depth=0.0,
            min_depth=0,
            ambiguous_bases=0,
            failure_reason=f"requires at least {min_depth} reads",
        )

    base_votes: list[Counter[str]] = [Counter() for _ in ref]
    base_observations = [0 for _ in ref]
    # Reads spanning a position, whatever the base quality. Quality decides
    # whether a base is *counted*, and must not decide how many reads were there:
    # shrinking the denominator raises the deletion share, and at nanopore
    # deletion rates that silently shortened 71 clones in a full run.
    spanning_reads = [0 for _ in ref]
    deletion_votes = [0 for _ in ref]
    insert_votes: dict[int, Counter[str]] = defaultdict(Counter)
    low_quality_bases = 0

    for read in selected:
        query = read.sequence.upper()
        if read.qualities is not None and len(read.qualities) != len(query):
            raise ValueError(f"quality length mismatch for {read.read_uid}")
        result = edlib.align(query, ref, mode="NW", task="path")
        cigar = result.get("cigar")
        if not cigar:
            raise RuntimeError(f"edlib produced no path for {read.read_uid}")
        qi = ri = 0
        for count_text, operation in _CIGAR_TOKEN.findall(cigar):
            count = int(count_text)
            if operation in "=X":
                for _ in range(count):
                    # Quality gates whether an observation is counted; it does not
                    # scale the vote. Weighting votes by Phred value made the
                    # support fraction a fraction of weighted votes, which is not
                    # the quantity a threshold on "support" appears to describe,
                    # and is not what a binomial test can be run on.
                    spanning_reads[ri] += 1
                    if (
                        read.qualities is None
                        or read.qualities[qi] >= min_base_quality
                    ):
                        base_votes[ri][query[qi]] += 1
                        base_observations[ri] += 1
                    else:
                        low_quality_bases += 1
                    qi += 1
                    ri += 1
            elif operation == "I":
                insertion = query[qi : qi + count]
                insert_votes[ri][insertion] += 1
                qi += count
            else:  # D: target/reference base absent from query
                for _ in range(count):
                    deletion_votes[ri] += 1
                    ri += 1

    # The group's own background error rate, from how often its reads disagree
    # with the reference at all. The mean rather than the median: the median sits
    # at the easy positions and understates what error does at the hard ones,
    # which made the test fire on systematic error.
    substitution_disagreement: list[float] = []
    deletion_disagreement: list[float] = []
    for ri, ref_base in enumerate(ref):
        observed = base_observations[ri]
        depth = observed + deletion_votes[ri]
        if observed:
            substitution_disagreement.append(
                (observed - base_votes[ri].get(ref_base, 0)) / observed
            )
        if depth:
            deletion_disagreement.append(deletion_votes[ri] / depth)

    def _rate(values: list[float]) -> float:
        mean = sum(values) / len(values) if values else _MIN_ERROR_RATE
        return min(max(mean, _MIN_ERROR_RATE), _MAX_ERROR_RATE)

    substitution_error = _rate(substitution_disagreement)
    deletion_error = _rate(deletion_disagreement)
    error_rate = substitution_error
    testable = sum(1 for value in substitution_disagreement if value > 0) or 1

    sequence_parts: list[str] = []
    depths: list[int] = []
    ambiguous = 0
    mixed: list[tuple[int, str, str, tuple[float, ...]]] = []
    weakest = 1.0
    for ri, ref_base in enumerate(ref):
        if ri in insert_votes:
            insertion, support = sorted(
                insert_votes[ri].items(), key=lambda item: (-item[1], item[0])
            )[0]
            insertion_depth = spanning_reads[ri] + deletion_votes[ri] + support
            if caller == "majority":
                accept = support >= min_depth and support / len(selected) >= min_support
            else:
                accept = (
                    support >= min_depth
                    and insertion_depth
                    and support / insertion_depth >= 0.5
                    and binomial_upper_tail(support, insertion_depth, deletion_error)
                    * testable
                    < significance
                )
            if accept:
                sequence_parts.append(insertion)
        votes = base_votes[ri]
        base_depth = base_observations[ri]
        position_depth = base_depth + deletion_votes[ri]
        depths.append(position_depth)
        if base_depth < min_depth and position_depth < min_depth:
            ambiguous += 1
            sequence_parts.append("N")
            continue
        alleles = Counter(votes)
        if deletion_votes[ri]:
            alleles["-"] = deletion_votes[ri]
        # Every read that reached this position, whatever its base quality.
        reads_spanning = spanning_reads[ri] + deletion_votes[ri]
        if caller == "majority":
            base, is_mixed, support = _majority_call(
                alleles, depth=position_depth, min_support=min_support
            )
            competing = ("", ())
        else:
            base, is_mixed, support, competing = _call_position(
                alleles,
                ref_base,
                depth=position_depth,
                spanning=reads_spanning,
                substitution_error=substitution_error,
                deletion_error=deletion_error,
                positions=testable,
                significance=significance,
                minor_fraction=minor_fraction,
                deletion_support=min_support,
            )
        if is_mixed:
            mixed.append((ri, ref_base, competing[0], competing[1]))
        if base == "N":
            ambiguous += 1
        weakest = min(weakest, support)
        if base != "-":
            sequence_parts.append(base)

    sequence = "".join(sequence_parts)
    status = "mixed_variants" if ambiguous else "consensus_pass"
    return ConsensusResult(
        sequence=sequence,
        status=status,
        n_reads_available=available,
        n_reads_used=len(selected),
        contributor_ids=tuple(r.read_uid for r in selected),
        mean_depth=sum(depths) / len(depths) if depths else math.nan,
        min_depth=min(depths) if depths else 0,
        ambiguous_bases=ambiguous,
        mixed_positions=len(mixed),
        mixed_alleles=tuple(mixed),
        background_error_rate=error_rate,
        deletion_error_rate=deletion_error,
        weakest_support=weakest,
        low_quality_bases=low_quality_bases,
    )


def _fasta_for_reads(reads: Sequence[ConsensusRead]) -> str:
    """Serialize selected reads with safe deterministic local identifiers."""

    parts: list[str] = []
    for index, read in enumerate(reads, start=1):
        sequence = read.sequence.upper()
        if not sequence:
            raise ValueError(f"consensus read {read.read_uid!r} has an empty sequence")
        if "\n" in sequence or "\r" in sequence:
            raise ValueError(
                f"consensus read {read.read_uid!r} contains a line break in its sequence"
            )
        if read.qualities is not None and len(read.qualities) != len(sequence):
            raise ValueError(f"quality length mismatch for {read.read_uid}")
        parts.append(f">read_{index:06d}\n{sequence}\n")
    return "".join(parts)


def _run_tool(
    name: str,
    command: Sequence[str],
    *,
    timeout_seconds: float,
) -> str:
    """Run one external tool without a shell and return its captured stdout."""

    try:
        completed = subprocess.run(
            list(command),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
            check=False,
            shell=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ConsensusBackendError(
            f"{name} timed out after {timeout_seconds:g} seconds while building consensus"
        ) from exc
    except OSError as exc:
        raise ConsensusBackendError(f"could not execute {name}: {exc}") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "no diagnostic output").strip()
        if len(detail) > 2_000:
            detail = detail[:2_000] + "…"
        raise ConsensusBackendError(
            f"{name} failed with exit code {completed.returncode}: {detail}"
        )
    if not completed.stdout.strip():
        raise ConsensusBackendError(f"{name} completed but produced no output")
    return completed.stdout


def _spoa_consensus_sequence(output: str) -> str:
    """Parse SPOA output exactly as the legacy notebook did, with validation."""

    sequence = "".join(
        line.strip() for line in output.splitlines() if not line.startswith(">")
    ).replace("-", "")
    sequence = sequence.upper()
    if not sequence:
        raise ConsensusBackendError("SPOA produced no consensus sequence")
    if set(sequence) - set("ACGTRYSWKMBDHVN"):
        invalid = "".join(sorted(set(sequence) - set("ACGTRYSWKMBDHVN")))
        raise ConsensusBackendError(
            f"SPOA consensus contains unsupported symbol(s): {invalid}"
        )
    return sequence


def _build_mafft_spoa_consensus(
    reference: str,
    reads: Sequence[ConsensusRead],
    *,
    group_id: str,
    min_depth: int,
    max_reads: int,
    seed: int,
    threads: int,
    mafft_path: str | Path | None,
    spoa_path: str | Path | None,
    timeout_seconds: float,
) -> ConsensusResult:
    """Reproduce the legacy MAFFT-to-SPOA command chain deterministically."""

    if not reference:
        raise ValueError("reference cannot be empty")
    if min_depth < 1 or max_reads < 1:
        raise ValueError("min_depth and max_reads must be positive")
    if threads < 1:
        raise ValueError("threads must be positive")
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be positive")

    available = len(reads)
    selected = select_reads(reads, max_reads=max_reads, seed=seed, group_id=group_id)
    contributor_ids = tuple(read.read_uid for read in selected)
    if len(selected) < min_depth:
        return ConsensusResult(
            sequence="",
            status="low_depth",
            n_reads_available=available,
            n_reads_used=len(selected),
            contributor_ids=contributor_ids,
            mean_depth=0.0,
            min_depth=0,
            ambiguous_bases=0,
            backend="mafft_spoa",
            failure_reason=f"requires at least {min_depth} reads",
        )

    mafft, spoa = require_mafft_spoa(mafft_path=mafft_path, spoa_path=spoa_path)
    with tempfile.TemporaryDirectory(prefix="nanopore3_mafft_spoa_") as temp_dir:
        directory = Path(temp_dir)
        reads_path = directory / "reads.fasta"
        alignment_path = directory / "alignment.fasta"
        reads_path.write_text(_fasta_for_reads(selected), encoding="ascii")

        alignment = _run_tool(
            "MAFFT",
            [mafft, "--quiet", "--thread", str(threads), str(reads_path)],
            timeout_seconds=timeout_seconds,
        )
        alignment_path.write_text(alignment, encoding="ascii")
        spoa_output = _run_tool(
            "SPOA",
            [spoa, "--algorithm", "msa", str(alignment_path)],
            timeout_seconds=timeout_seconds,
        )

    sequence = _spoa_consensus_sequence(spoa_output)
    ambiguous = sequence.count("N")
    return ConsensusResult(
        sequence=sequence,
        status="mixed_variants" if ambiguous else "consensus_pass",
        n_reads_available=available,
        n_reads_used=len(selected),
        contributor_ids=contributor_ids,
        mean_depth=math.nan,
        min_depth=0,
        ambiguous_bases=ambiguous,
        backend="mafft_spoa",
    )


def build_reference_consensus(
    reference: str,
    reads: Sequence[ConsensusRead],
    *,
    group_id: str,
    min_depth: int = 3,
    max_reads: int = 100,
    min_support: float = 0.60,
    seed: int = 0,
    min_base_quality: int = 10,
    significance: float = 0.05,
    minor_fraction: float = 0.20,
    caller: str = "statistical",
    backend: str = "portable",
    threads: int = 1,
    mafft_path: str | Path | None = None,
    spoa_path: str | Path | None = None,
    timeout_seconds: float = 3_600,
) -> ConsensusResult:
    """Build a consensus with an explicitly selected deterministic backend.

    ``portable`` preserves the reference-guided edlib pileup. ``mafft_spoa``
    reproduces the legacy external command chain, while using stable hash-based
    contributor selection rather than FASTQ order.
    """

    selected_backend = validate_consensus_backend(backend)
    if selected_backend == "portable":
        return _build_portable_consensus(
            reference,
            reads,
            group_id=group_id,
            min_depth=min_depth,
            max_reads=max_reads,
            min_support=min_support,
            seed=seed,
            min_base_quality=min_base_quality,
            significance=significance,
            minor_fraction=minor_fraction,
            caller=caller,
        )
    return _build_mafft_spoa_consensus(
        reference,
        reads,
        group_id=group_id,
        min_depth=min_depth,
        max_reads=max_reads,
        seed=seed,
        threads=threads,
        mafft_path=mafft_path,
        spoa_path=spoa_path,
        timeout_seconds=timeout_seconds,
    )
