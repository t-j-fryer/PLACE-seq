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

    def to_dict(self) -> dict[str, object]:
        value = asdict(self)
        value.update({f"alignment_{k}": v for k, v in asdict(self.metrics).items()})
        del value["metrics"]
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
) -> QcResult:
    """Evaluate independent criteria without treating missing boundaries as failure.

    Supplying the constant sequence that flanks the consensus in the full open
    reading frame enables the coding checks.  ``upstream_constant`` must begin at
    the start codon; the reading frame is then defined by construction rather
    than assumed, and the assembled
    ``upstream_constant + consensus + downstream_constant`` is translated.
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
    if upstream_constant is not None and downstream_constant is not None and consensus:
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
    )

