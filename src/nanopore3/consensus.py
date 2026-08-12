"""Deterministic consensus construction with optional external-tool backends.

The default ``portable`` backend deliberately uses only edlib.  ``mafft_spoa`` is
an opt-in compatibility backend for the legacy Nanopore2 workflow; importing and
using the package never requires either external executable unless that backend is
selected explicitly.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Iterable, Sequence

import edlib


_CIGAR_TOKEN = re.compile(r"(\d+)([=XID])")
SUPPORTED_CONSENSUS_BACKENDS = frozenset({"portable", "mafft_spoa"})


class ConsensusBackendUnavailable(RuntimeError):
    """Raised when an explicitly selected optional backend cannot be executed."""


class ConsensusBackendError(RuntimeError):
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


def _stable_rank(read_uid: str, seed: int, group_id: str) -> str:
    value = f"{seed}\0{group_id}\0{read_uid}".encode("utf-8")
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


def _weight(qualities: tuple[int, ...] | None, index: int) -> int:
    if qualities is None:
        return 1
    return max(1, qualities[index])


def _choice(votes: Counter[str], *, min_support: float) -> str:
    if not votes:
        return "N"
    ordered = sorted(votes.items(), key=lambda item: (-item[1], item[0]))
    winner, support = ordered[0]
    total = sum(votes.values())
    if total == 0 or support / total < min_support:
        return "N"
    return winner


def _build_portable_consensus(
    reference: str,
    reads: Sequence[ConsensusRead],
    *,
    group_id: str,
    min_depth: int = 3,
    max_reads: int = 100,
    min_support: float = 0.60,
    seed: int = 0,
) -> ConsensusResult:
    """Build the edlib reference-guided consensus used by the portable backend."""

    ref = reference.upper()
    if not ref:
        raise ValueError("reference cannot be empty")
    if min_depth < 1 or max_reads < 1:
        raise ValueError("min_depth and max_reads must be positive")
    if not 0.5 <= min_support <= 1.0:
        raise ValueError("min_support must be between 0.5 and 1.0")

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
    deletion_votes = [0 for _ in ref]
    insert_votes: dict[int, Counter[str]] = defaultdict(Counter)

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
                    base_votes[ri][query[qi]] += _weight(read.qualities, qi)
                    base_observations[ri] += 1
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

    sequence_parts: list[str] = []
    depths: list[int] = []
    ambiguous = 0
    for ri, _ref_base in enumerate(ref):
        if ri in insert_votes:
            insertion, support = sorted(
                insert_votes[ri].items(), key=lambda item: (-item[1], item[0])
            )[0]
            if support >= min_depth and support / len(selected) >= min_support:
                sequence_parts.append(insertion)
        votes = base_votes[ri]
        base_depth = base_observations[ri]
        position_depth = base_depth + deletion_votes[ri]
        depths.append(position_depth)
        if position_depth and deletion_votes[ri] / position_depth >= min_support:
            continue
        base = "N" if base_depth < min_depth else _choice(votes, min_support=min_support)
        if base == "N":
            ambiguous += 1
        sequence_parts.append(base)

    sequence = "".join(sequence_parts)
    status = "heterogeneous" if ambiguous else "consensus_pass"
    return ConsensusResult(
        sequence=sequence,
        status=status,
        n_reads_available=available,
        n_reads_used=len(selected),
        contributor_ids=tuple(r.read_uid for r in selected),
        mean_depth=sum(depths) / len(depths) if depths else math.nan,
        min_depth=min(depths) if depths else 0,
        ambiguous_bases=ambiguous,
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
        status="heterogeneous" if ambiguous else "consensus_pass",
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
