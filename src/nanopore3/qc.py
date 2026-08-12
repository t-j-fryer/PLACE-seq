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

    def to_dict(self) -> dict[str, object]:
        value = asdict(self)
        value.update({f"alignment_{k}": v for k, v in asdict(self.metrics).items()})
        del value["metrics"]
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


def _translate_has_internal_stop(sequence: str, frame: int = 0) -> bool:
    stops = {"TAA", "TAG", "TGA"}
    codons = [sequence[i : i + 3] for i in range(frame, len(sequence) - 2, 3)]
    return any(codon in stops for codon in codons[:-1])


def evaluate_consensus(
    consensus: str,
    reference: str,
    *,
    min_identity: float = 0.98,
    min_query_coverage: float = 0.95,
    min_reference_coverage: float = 0.95,
    length_tolerance: int = 10,
    coding_start: int | None = None,
    coding_end: int | None = None,
) -> QcResult:
    """Evaluate independent criteria without treating missing boundaries as failure."""

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

    if coding_start is None or coding_end is None:
        frame: QcState = "not_evaluable"
        stops: QcState = "not_evaluable"
    elif not (0 <= coding_start < coding_end <= len(consensus)):
        frame = stops = "not_evaluable"
    else:
        coding = consensus[coding_start:coding_end].upper()
        frame = "pass" if len(coding) % 3 == 0 else "fail"
        stops = "fail" if _translate_has_internal_stop(coding) else "pass"

    states = (full, expected, frame, stops)
    evaluable = tuple(state for state in states if state != "not_evaluable")
    overall: QcState = "fail" if "fail" in evaluable else ("pass" if evaluable else "not_evaluable")
    reason = "one or more criteria failed" if overall == "fail" else "all evaluable criteria passed"
    return QcResult(full, expected, frame, stops, overall, metrics, reason)

