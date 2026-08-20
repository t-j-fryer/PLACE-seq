"""Small, explicit QC primitives with pass/fail/not-evaluable semantics."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal

import edlib

QcState = Literal["pass", "fail", "not_evaluable"]


@dataclass(frozen=True, slots=True)
class AlignmentMetrics:
    edit_distance: int
    identity: float
    query_coverage: float
    reference_coverage: float


@dataclass(frozen=True, slots=True)
class RegionMetrics:
    """Edit distance and identity over one named span of the reference."""

    name: str
    reference_length: int
    edit_distance: int
    identity: float


@dataclass(frozen=True, slots=True)
class QcResult:
    full_amplicon: QcState
    expected_length: QcState
    reading_frame: QcState
    internal_stops: QcState
    overall: QcState
    metrics: AlignmentMetrics
    reason: str
    # Residues before the terminal stop, and the 0-based codon index of the
    # first internal stop. Both are None when the coding checks did not run.
    protein_length: int | None = None
    internal_stop_codon: int | None = None
    # Per-region accuracy, present only when the reference declares spans.
    regions: tuple[RegionMetrics, ...] = ()

    def to_dict(self) -> dict[str, object]:
        value = asdict(self)
        value.update({f"alignment_{k}": v for k, v in asdict(self.metrics).items()})
        del value["metrics"]
        for region in self.regions:
            value[f"{region.name}_edit_distance"] = region.edit_distance
            value[f"{region.name}_identity"] = f"{region.identity:.6f}"
            value[f"{region.name}_length"] = region.reference_length
        del value["regions"]
        for key in ("protein_length", "internal_stop_codon"):
            if value[key] is None:
                value[key] = ""
        return value


def global_alignment_metrics(query: str, reference: str) -> AlignmentMetrics:
    """Return intentionally simple whole-sequence edit and coverage metrics."""

    q, r = query.upper(), reference.upper()
    if not q or not r:
        return AlignmentMetrics(-1, 0.0, 0.0, 0.0)
    distance = edlib.align(q, r, mode="NW", task="distance")["editDistance"]
    denominator = max(len(q), len(r))
    return AlignmentMetrics(
        edit_distance=distance,
        identity=max(0.0, 1.0 - distance / denominator),
        query_coverage=min(1.0, len(r) / len(q)),
        reference_coverage=min(1.0, len(q) / len(r)),
    )


STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

# Standard genetic code. Any codon containing an ambiguity code translates to
# "X" rather than guessing, so a low-support consensus base cannot silently
# become a confident amino acid.
_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(sequence: str) -> str:
    """Translate frame 0, mapping any codon with an ambiguity code to ``X``."""

    upper = sequence.upper()
    return "".join(
        _CODON_TABLE.get(upper[index : index + 3], "X")
        for index in range(0, len(upper) - len(upper) % 3, 3)
    )


def first_internal_stop(protein: str) -> int | None:
    """Return the codon index of the first stop before the final codon."""

    for index, residue in enumerate(protein[:-1] if protein else protein):
        if residue == "*":
            return index
    return None


def _cigar_steps(cigar: str) -> list[tuple[int, str]]:
    steps: list[tuple[int, str]] = []
    number = ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        steps.append((int(number or 1), char))
        number = ""
    return steps


def region_metrics(
    consensus: str,
    reference: str,
    spans: dict[str, tuple[int, int]],
) -> tuple[RegionMetrics, ...]:
    """Edit distance per named span of the reference.

    A full-length amplicon is not uniform: the insert is the designed part and
    the flanks are the vector, and an error in one means something different from
    an error in the other.  A single whole-sequence identity hides that, so the
    alignment is walked once and each edit charged to the reference span it falls
    in.  Insertions are charged to the span at the current reference position.
    """

    query, target = consensus.upper(), reference.upper()
    if not query or not target or not spans:
        return ()
    result = edlib.align(query, target, mode="NW", task="path")
    cigar = result.get("cigar") or ""
    edits: dict[str, int] = {name: 0 for name in spans}
    position = 0  # reference coordinate
    for count, operation in _cigar_steps(cigar):
        for _ in range(count):
            if operation in ("X", "I", "D"):
                for name, (start, end) in spans.items():
                    if start <= min(position, len(target) - 1) < end:
                        edits[name] += 1
                        break
            if operation != "I":  # "=", "X" and "D" consume the reference
                position += 1
    return tuple(
        RegionMetrics(
            name=name,
            reference_length=end - start,
            edit_distance=edits[name],
            identity=max(0.0, 1.0 - edits[name] / (end - start)) if end > start else 0.0,
        )
        for name, (start, end) in spans.items()
    )


def _locate(motif: str, sequence: str, max_edits: int) -> tuple[int, int] | None:
    """Best fuzzy location of a motif, as half-open coordinates."""

    if not motif or not sequence:
        return None
    result = edlib.align(motif.upper(), sequence.upper(), mode="HW", task="locations")
    if result["editDistance"] < 0 or result["editDistance"] > max_edits:
        return None
    start, end = result["locations"][0]
    return start, end + 1


def coding_checks_in_consensus(
    consensus: str,
    *,
    start_anchor: str,
    stop_anchor: str,
    max_anchor_edits: int = 3,
) -> tuple[QcState, QcState, str, int | None, str]:
    """Frame and internal-stop checks for a consensus that spans the whole ORF.

    In insert mode the ORF is assembled by prepending and appending the constant
    regions.  Here the consensus already contains them, so the ORF is located
    instead: from the start anchor - which must begin at the start codon - to the
    end of the stop anchor, which ends at the terminal stop.  Frame is then the
    measured distance between them, not an assumption.
    """

    sequence = consensus.upper()
    start = _locate(start_anchor, sequence, max_anchor_edits)
    stop = _locate(stop_anchor, sequence, max_anchor_edits)
    if start is None or stop is None:
        missing = "start codon" if start is None else "terminal stop"
        return "not_evaluable", "not_evaluable", "", None, f"{missing} anchor not found"
    orf = sequence[start[0] : stop[1]]
    remainder = len(orf) % 3
    if remainder:
        return (
            "fail",
            "not_evaluable",
            "",
            None,
            f"open reading frame is {len(orf)} nt, {remainder} past a codon boundary",
        )
    protein = translate(orf)
    index = first_internal_stop(protein)
    detail = f"internal stop codon at residue {index + 1}" if index is not None else ""
    return "pass", ("fail" if index is not None else "pass"), protein, index, detail


def evaluate_consensus(
    consensus: str,
    reference: str,
    *,
    min_identity: float = 0.98,
    min_query_coverage: float = 0.95,
    min_reference_coverage: float = 0.95,
    length_tolerance: int = 10,
    upstream_constant: str | None = None,
    downstream_constant: str | None = None,
    coding_anchors: tuple[str, str] | None = None,
    spans: dict[str, tuple[int, int]] | None = None,
) -> QcResult:
    """Evaluate independent criteria without treating missing boundaries as failure.

    Two ways to enable the coding checks, for the two shapes of consensus:

    * insert mode - supply ``upstream_constant`` and ``downstream_constant``, the
      sequence that flanks the consensus in the full open reading frame.  The
      upstream constant must begin at the start codon, and
      ``upstream + consensus + downstream`` is translated.
    * full-length mode - supply ``coding_anchors`` as (start, stop), where the
      consensus already contains the whole ORF.  The anchors are located in the
      consensus, so frame is measured rather than assumed.

    ``spans`` names regions of the reference to report accuracy over separately,
    which is what distinguishes an error in a designed insert from one in the
    vector around it.
    """

    metrics = global_alignment_metrics(consensus, reference)
    full: QcState = (
        "pass"
        if metrics.identity >= min_identity
        and metrics.query_coverage >= min_query_coverage
        and metrics.reference_coverage >= min_reference_coverage
        else "fail"
    )
    expected: QcState = (
        "pass" if abs(len(consensus) - len(reference)) <= length_tolerance else "fail"
    )

    frame: QcState = "not_evaluable"
    stops: QcState = "not_evaluable"
    protein = ""
    stop_codon_index: int | None = None
    detail = ""
    if coding_anchors is not None and consensus:
        # Full-length mode: the ORF is inside the consensus, so locate it rather
        # than assembling one around it.
        frame, stops, protein, stop_codon_index, detail = coding_checks_in_consensus(
            consensus, start_anchor=coding_anchors[0], stop_anchor=coding_anchors[1]
        )
    elif upstream_constant is not None and downstream_constant is not None and consensus:
        orf = (upstream_constant + consensus + downstream_constant).upper()
        remainder = len(orf) % 3
        frame = "pass" if remainder == 0 else "fail"
        if remainder:
            # A frameshifted ORF translates to noise, so reporting stop codons
            # from it would be misleading rather than informative.
            detail = (
                f"assembled reading frame is {len(orf)} nt, "
                f"{remainder} past a codon boundary"
            )
        else:
            protein = translate(orf)
            stop_codon_index = first_internal_stop(protein)
            stops = "fail" if stop_codon_index is not None else "pass"
            if stop_codon_index is not None:
                detail = f"internal stop codon at residue {stop_codon_index + 1}"

    states = (full, expected, frame, stops)
    evaluable = tuple(state for state in states if state != "not_evaluable")
    overall: QcState = "fail" if "fail" in evaluable else ("pass" if evaluable else "not_evaluable")
    if overall == "fail":
        reason = detail or "one or more criteria failed"
    else:
        reason = "all evaluable criteria passed"
    return QcResult(
        full, expected, frame, stops, overall, metrics, reason,
        protein_length=len(protein.rstrip("*")) if protein else None,
        internal_stop_codon=stop_codon_index,
        regions=region_metrics(consensus, reference, spans) if spans else (),
    )

