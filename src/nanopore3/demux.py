"""Evidence-preserving barcode validation and calling."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

import edlib

from .errors import Nanopore3Error
from .sequence import edlib_iupac_equalities, normalize_sequence, reverse_complement

BarcodeStatus = Literal["assigned", "ambiguous", "unassigned", "conflicting"]
ReadOrientation = Literal["forward", "reverse", "unknown"]
ReadEnd = Literal["head", "tail"]


class BarcodeValidationError(Nanopore3Error, ValueError):
    """Raised when a barcode panel is invalid or not safely distinguishable."""


@dataclass(frozen=True, slots=True)
class BarcodeEvidence:
    """The best semi-global match for one barcode form in one read-end window."""

    barcode_id: str
    distance: int
    read_end: ReadEnd
    orientation: Literal["forward", "reverse"]
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class BarcodeCall:
    """A barcode decision with sufficient evidence for later audit."""

    status: BarcodeStatus
    barcode_id: str | None
    orientation: ReadOrientation
    best_distance: int | None
    second_distance: int | None
    margin: int | None
    matched_end: ReadEnd | None
    start: int | None
    end: int | None
    reason: str
    evidence: tuple[BarcodeEvidence, ...]


@dataclass(frozen=True, slots=True)
class PreparedBarcodePanel:
    """Normalized immutable barcode forms safe to reuse across all reads."""

    sequences: tuple[tuple[str, str, str], ...]
    uses_iupac: bool


def prepare_barcode_panel(barcodes: Mapping[str, str]) -> PreparedBarcodePanel:
    """Normalize a panel and precompute reverse complements exactly once."""

    normalized = _normalized_barcodes(barcodes)
    sequences = tuple(
        (barcode_id, sequence, reverse_complement(sequence))
        for barcode_id, sequence in normalized.items()
    )
    return PreparedBarcodePanel(
        sequences=sequences,
        uses_iupac=any(set(sequence) - set("ACGT") for sequence in normalized.values()),
    )


def _normalized_barcodes(barcodes: Mapping[str, str]) -> dict[str, str]:
    if not barcodes:
        raise BarcodeValidationError("at least one barcode is required")
    result: dict[str, str] = {}
    for raw_id, raw_sequence in barcodes.items():
        barcode_id = str(raw_id).strip()
        if not barcode_id:
            raise BarcodeValidationError("barcode identifiers must not be empty")
        if barcode_id in result:
            raise BarcodeValidationError(f"duplicate barcode identifier: {barcode_id!r}")
        result[barcode_id] = normalize_sequence(raw_sequence)
    return dict(sorted(result.items()))


def _global_distance(left: str, right: str) -> int:
    result = edlib.align(
        left,
        right,
        mode="NW",
        task="distance",
        additionalEqualities=edlib_iupac_equalities(),
    )
    return int(result["editDistance"])


def validate_barcodes(
    barcodes: Mapping[str, str],
    *,
    max_edits: int,
    min_margin: int = 1,
    orientation_aware: bool = True,
) -> dict[str, str]:
    """Normalize a barcode panel and reject intrinsically unsafe pairs.

    A pair must be separated by more than ``2 * max_edits + min_margin - 1``
    edits, which prevents their accepted error neighbourhoods from violating the
    requested score margin.  Reverse-complement forms are included when reads
    will be searched in both orientations.
    """

    if max_edits < 0:
        raise ValueError("max_edits must be non-negative")
    if min_margin < 0:
        raise ValueError("min_margin must be non-negative")
    normalized = _normalized_barcodes(barcodes)
    minimum_safe_distance = 2 * max_edits + min_margin
    identifiers = tuple(normalized)
    unsafe: list[str] = []
    for index, left_id in enumerate(identifiers):
        for right_id in identifiers[index + 1 :]:
            left = normalized[left_id]
            right = normalized[right_id]
            distances = [_global_distance(left, right)]
            if orientation_aware:
                distances.extend(
                    (
                        _global_distance(left, reverse_complement(right)),
                        _global_distance(reverse_complement(left), right),
                    )
                )
            distance = min(distances)
            if distance < minimum_safe_distance:
                unsafe.append(f"{left_id!r}/{right_id!r} (distance {distance})")
    if unsafe:
        pairs = "; ".join(unsafe)
        raise BarcodeValidationError(
            f"barcode panel is not distinguishable at max_edits={max_edits}, "
            f"min_margin={min_margin}: {pairs}"
        )
    return normalized


def _best_location(
    pattern: str,
    target: str,
    *,
    uses_iupac: bool,
    target_uses_iupac: bool,
    max_edits: int | None = None,
) -> tuple[int, int, int] | None:
    kwargs = {}
    if uses_iupac or target_uses_iupac:
        kwargs["additionalEqualities"] = edlib_iupac_equalities()
    if max_edits is not None:
        kwargs["k"] = max_edits
    result = edlib.align(
        pattern,
        target,
        mode="HW",
        task="locations",
        **kwargs,
    )
    distance = int(result["editDistance"])
    locations = [
        location
        for location in result.get("locations", ())
        if location[0] is not None
        and location[1] is not None
        and int(location[0]) >= 0
        and int(location[1]) >= int(location[0])
    ]
    if distance < 0 or not locations:
        return None
    start, end = min((int(start), int(end)) for start, end in locations)
    return distance, start, end


def _end_evidence(
    read: str,
    panel: PreparedBarcodePanel,
    window_size: int,
    search_ends: tuple[ReadEnd, ...],
    allow_reverse_complement: bool,
    max_edits: int | None,
) -> tuple[BarcodeEvidence, ...]:
    all_windows: tuple[tuple[ReadEnd, str, int], ...] = (
        ("head", read[:window_size], 0),
        ("tail", read[-window_size:], max(0, len(read) - window_size)),
    )
    windows = tuple(
        (read_end, window, offset, bool(set(window) - set("ACGT")))
        for read_end, window, offset in all_windows
        if read_end in search_ends
    )
    evidence: list[BarcodeEvidence] = []
    for barcode_id, barcode, reverse in panel.sequences:
        for read_end, window, offset, window_uses_iupac in windows:
            # A forward barcode at the head or its reverse complement at the
            # tail supports forward read orientation; the inverse supports reverse.
            forms: tuple[tuple[str, Literal["forward", "reverse"]], ...]
            if read_end == "head":
                forms = ((barcode, "forward"), (reverse, "reverse"))
            else:
                forms = ((reverse, "forward"), (barcode, "reverse"))
            if not allow_reverse_complement:
                forms = tuple(form for form in forms if form[1] == "forward")
            for pattern, orientation in forms:
                match = _best_location(
                    pattern,
                    window,
                    uses_iupac=panel.uses_iupac,
                    target_uses_iupac=window_uses_iupac,
                    max_edits=max_edits,
                )
                if match is None:
                    continue
                distance, start, end = match
                evidence.append(
                    BarcodeEvidence(
                        barcode_id=barcode_id,
                        distance=distance,
                        read_end=read_end,
                        orientation=orientation,
                        start=start + offset,
                        end=end + offset,
                    )
                )
    return tuple(
        sorted(
            evidence,
            key=lambda item: (
                item.distance,
                item.barcode_id,
                item.orientation,
                item.read_end,
                item.start,
                item.end,
            ),
        )
    )


def _representatives(evidence: Sequence[BarcodeEvidence]) -> tuple[BarcodeEvidence, ...]:
    """Retain one deterministic best hit per barcode, end, and orientation."""

    best: dict[tuple[str, ReadEnd, str], BarcodeEvidence] = {}
    for hit in evidence:
        key = (hit.barcode_id, hit.read_end, hit.orientation)
        if key not in best:
            best[key] = hit
    return tuple(
        sorted(
            best.values(),
            key=lambda item: (
                item.distance,
                item.barcode_id,
                item.orientation,
                item.read_end,
                item.start,
            ),
        )
    )


def call_barcode(
    sequence: str,
    barcodes: Mapping[str, str],
    *,
    window_size: int = 400,
    max_edits: int = 4,
    min_margin: int = 1,
    search_ends: tuple[ReadEnd, ...] = ("head", "tail"),
    allow_reverse_complement: bool = True,
    decision_policy: str = "best_margin",
    prepared_panel: PreparedBarcodePanel | None = None,
) -> BarcodeCall:
    """Call a barcode from both read ends using IUPAC-aware semi-global matches.

    The function does not assume that the barcode panel has passed strict
    preflight, allowing close panels to yield explicit ``ambiguous`` calls.
    Production callers should invoke :func:`validate_barcodes` once beforehand.
    """

    if window_size < 1:
        raise ValueError("window_size must be at least 1")
    if max_edits < 0 or min_margin < 0:
        raise ValueError("max_edits and min_margin must be non-negative")
    if decision_policy not in {
        "best_margin",
        "legacy_unique_threshold",
        "legacy_unique_best",
    }:
        raise ValueError(
            "decision_policy must be best_margin, legacy_unique_threshold, "
            "or legacy_unique_best"
        )
    if not search_ends or set(search_ends) - {"head", "tail"}:
        raise ValueError("search_ends must contain only 'head' and/or 'tail'")
    read = normalize_sequence(sequence)
    panel = prepared_panel or prepare_barcode_panel(barcodes)
    all_evidence = _representatives(
        _end_evidence(
            read,
            panel,
            min(window_size, len(read)),
            search_ends,
            allow_reverse_complement,
            max_edits if decision_policy.startswith("legacy_") else None,
        )
    )
    eligible = tuple(hit for hit in all_evidence if hit.distance <= max_edits)
    if not eligible:
        best_distance = all_evidence[0].distance if all_evidence else None
        return BarcodeCall(
            "unassigned", None, "unknown", best_distance, None, None, None, None, None,
            "no barcode met the maximum edit distance", all_evidence,
        )

    if decision_policy.startswith("legacy_"):
        owners = sorted({hit.barcode_id for hit in eligible})
        best = eligible[0]
        alternatives = [hit for hit in all_evidence if hit.barcode_id != best.barcode_id]
        second = alternatives[0] if alternatives else None
        second_distance = None if second is None else second.distance
        margin = (
            None if second_distance is None else second_distance - best.distance
        )
        best_owners = sorted(
            {hit.barcode_id for hit in eligible if hit.distance == best.distance}
        )
        ambiguous = (
            len(owners) != 1
            if decision_policy == "legacy_unique_threshold"
            else len(best_owners) != 1
        )
        if ambiguous:
            return BarcodeCall(
                "ambiguous",
                None,
                "unknown",
                best.distance,
                second_distance,
                margin,
                best.read_end,
                best.start,
                best.end,
                (
                    "multiple barcode identities met the legacy edit threshold"
                    if decision_policy == "legacy_unique_threshold"
                    else "multiple barcode identities tied at the best legacy distance"
                ),
                all_evidence,
            )
        owner = owners[0] if decision_policy == "legacy_unique_threshold" else best_owners[0]
        owner_hits = [hit for hit in eligible if hit.barcode_id == owner]
        # The legacy plate code resolves equal head/tail evidence in favour of
        # the tail (``tail_best <= head_best``), retaining the original read.
        owner_best = min(
            owner_hits,
            key=lambda hit: (
                hit.distance,
                0 if hit.read_end == "tail" else 1,
                hit.start,
                hit.orientation,
            ),
        )
        # Nanopore2 plate demultiplexing defined orientation by the winning end:
        # tail evidence retained the read, while head evidence reverse-complemented it.
        orientation = "forward" if owner_best.read_end == "tail" else "reverse"
        if not allow_reverse_complement:
            orientation = "forward"
        return BarcodeCall(
            "assigned",
            owner,
            orientation,
            owner_best.distance,
            second_distance,
            margin,
            owner_best.read_end,
            owner_best.start,
            owner_best.end,
            (
                "exactly one barcode identity met the legacy edit threshold"
                if decision_policy == "legacy_unique_threshold"
                else "one barcode identity had the unique best legacy distance"
            ),
            all_evidence,
        )

    # Each end may independently establish a confident identity.  Conflicting
    # identities or orientations are reported before pooling all evidence.
    confident_end_hits: list[BarcodeEvidence] = []
    for read_end in ("head", "tail"):
        end_hits = sorted(
            (hit for hit in eligible if hit.read_end == read_end),
            key=lambda item: (item.distance, item.barcode_id, item.orientation, item.start),
        )
        if not end_hits:
            continue
        best = end_hits[0]
        alternatives = [
            hit
            for hit in end_hits[1:]
            if (hit.barcode_id, hit.orientation) != (best.barcode_id, best.orientation)
        ]
        second = alternatives[0] if alternatives else None
        if second is None or second.distance - best.distance >= min_margin:
            confident_end_hits.append(best)

    if len(confident_end_hits) == 2:
        head, tail = confident_end_hits
        if (head.barcode_id, head.orientation) != (tail.barcode_id, tail.orientation):
            return BarcodeCall(
                "conflicting", None, "unknown", min(head.distance, tail.distance),
                max(head.distance, tail.distance), abs(tail.distance - head.distance), None,
                None, None, "read ends support conflicting barcode identities or orientations",
                all_evidence,
            )

    best = eligible[0]
    alternatives = [
        hit
        for hit in eligible[1:]
        if (hit.barcode_id, hit.orientation) != (best.barcode_id, best.orientation)
    ]
    second = alternatives[0] if alternatives else None
    second_distance = second.distance if second else None
    margin = second_distance - best.distance if second_distance is not None else None
    if second is not None and margin is not None and margin < min_margin:
        return BarcodeCall(
            "ambiguous", None, "unknown", best.distance, second_distance, margin, best.read_end,
            best.start, best.end, "best barcode/orientation does not meet the score margin",
            all_evidence,
        )
    return BarcodeCall(
        "assigned", best.barcode_id, best.orientation, best.distance, second_distance, margin,
        best.read_end, best.start, best.end, "barcode met distance and margin thresholds",
        all_evidence,
    )
