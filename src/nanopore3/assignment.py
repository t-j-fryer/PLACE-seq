"""Motif extraction and conservative reference assignment primitives."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Literal, Mapping

import edlib

from .sequence import (
    canonical_unique_kmers,
    edlib_iupac_equalities,
    normalize_sequence,
    reverse_complement,
)

AssignmentStatus = Literal[
    "assigned_unique",
    "assigned_alias_set",
    "ambiguous",
    "no_match",
    "motif_missing",
]


@dataclass(frozen=True, slots=True)
class MotifMatch:
    motif: str
    distance: int
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class InsertExtraction:
    status: Literal["found", "ambiguous", "motif_missing"]
    sequence: str | None
    orientation: Literal["forward", "reverse", "unknown"]
    left: MotifMatch | None
    right: MotifMatch | None
    reason: str


@dataclass(frozen=True, slots=True)
class CandidateScore:
    aliases: tuple[str, ...]
    kmer_score: float
    matching_kmers: int


@dataclass(frozen=True, slots=True)
class AlignmentEvidence:
    aliases: tuple[str, ...]
    kmer_score: float
    matching_kmers: int
    edit_distance: int
    identity: float
    query_coverage: float
    reference_coverage: float
    query_start: int
    query_end: int
    reference_start: int
    reference_end: int


@dataclass(frozen=True, slots=True)
class AssignmentCall:
    status: AssignmentStatus
    reference_ids: tuple[str, ...]
    orientation: Literal["forward", "reverse", "unknown"]
    query_sequence: str | None
    best: AlignmentEvidence | None
    second: AlignmentEvidence | None
    identity_margin: float | None
    reason: str
    candidates: tuple[CandidateScore, ...]
    alignments: tuple[AlignmentEvidence, ...]
    extraction: InsertExtraction | None


def _motif_locations(sequence: str, motif: str, max_edits: int) -> tuple[MotifMatch, ...]:
    result = edlib.align(
        motif,
        sequence,
        mode="HW",
        task="locations",
        k=max_edits,
        additionalEqualities=edlib_iupac_equalities(),
    )
    distance = int(result["editDistance"])
    if distance < 0:
        return ()
    matches = {
        (int(start), int(end))
        for start, end in result.get("locations", ())
        if start is not None and end is not None
    }
    return tuple(
        MotifMatch(motif, distance, start, end)
        for start, end in sorted(matches)
    )


def _extract_one_orientation(
    sequence: str,
    left_motif: str,
    right_motif: str,
    max_edits: int,
    min_insert_length: int,
    max_insert_length: int | None,
    expected_insert_length: int | None,
) -> tuple[tuple[int, int, int, int, MotifMatch, MotifMatch], ...]:
    left_matches = _motif_locations(sequence, left_motif, max_edits)
    right_matches = _motif_locations(sequence, right_motif, max_edits)
    pairs: list[tuple[int, int, int, int, MotifMatch, MotifMatch]] = []
    for left in left_matches:
        insert_start = left.end + 1
        for right in right_matches:
            insert_end = right.start
            length = insert_end - insert_start
            if length < min_insert_length:
                continue
            if max_insert_length is not None and length > max_insert_length:
                continue
            length_deviation = (
                abs(length - expected_insert_length) if expected_insert_length is not None else 0
            )
            # The final coordinate terms ensure deterministic choice among repeats.
            pairs.append(
                (
                    left.distance + right.distance,
                    length_deviation,
                    left.start,
                    right.start,
                    left,
                    right,
                )
            )
    return tuple(sorted(pairs, key=lambda pair: pair[:4]))


def extract_insert(
    sequence: str,
    left_motif: str,
    right_motif: str,
    *,
    max_edits: int = 2,
    search_reverse_complement: bool = True,
    min_insert_length: int = 0,
    max_insert_length: int | None = None,
    expected_insert_length: int | None = None,
) -> InsertExtraction:
    """Extract sequence between an ordered fuzzy motif pair.

    The best pair minimizes total motif edits, then expected-length deviation,
    then coordinates. Equal-scoring forward/reverse solutions with different
    inserts are rejected as ambiguous rather than chosen by processing order.
    """

    if max_edits < 0 or min_insert_length < 0:
        raise ValueError("edit and length thresholds must be non-negative")
    if max_insert_length is not None and max_insert_length < min_insert_length:
        raise ValueError("max_insert_length must be at least min_insert_length")
    read = normalize_sequence(sequence)
    left = normalize_sequence(left_motif)
    right = normalize_sequence(right_motif)
    orientations: list[tuple[Literal["forward", "reverse"], str]] = [("forward", read)]
    if search_reverse_complement:
        orientations.append(("reverse", reverse_complement(read)))

    solutions: list[
        tuple[
            tuple[int, int, int, int],
            Literal["forward", "reverse"],
            str,
            MotifMatch,
            MotifMatch,
        ]
    ] = []
    for orientation, oriented_read in orientations:
        pairs = _extract_one_orientation(
            oriented_read,
            left,
            right,
            max_edits,
            min_insert_length,
            max_insert_length,
            expected_insert_length,
        )
        for total_edits, length_deviation, left_start, right_start, left_hit, right_hit in pairs:
            key = (total_edits, length_deviation, left_start, right_start)
            insert = oriented_read[left_hit.end + 1 : right_hit.start]
            solutions.append((key, orientation, insert, left_hit, right_hit))

    if not solutions:
        return InsertExtraction(
            "motif_missing", None, "unknown", None, None,
            "no ordered motif pair met the edit and insert-length constraints",
        )
    solutions.sort(key=lambda item: (item[0], item[1], item[2]))
    best_key = solutions[0][0][:2]
    equally_scored = [solution for solution in solutions if solution[0][:2] == best_key]
    distinct = {(solution[1], solution[2]) for solution in equally_scored}
    if len(distinct) > 1:
        return InsertExtraction(
            "ambiguous", None, "unknown", None, None,
            "multiple orientations or motif pairs have equal evidence",
        )
    _, orientation, insert, left_hit, right_hit = equally_scored[0]
    return InsertExtraction(
        "found", insert, orientation, left_hit, right_hit,
        "ordered motif pair found",
    )


class ReferenceIndex:
    """Immutable reference alias groups and specificity-weighted k-mer index."""

    def __init__(
        self,
        references: Mapping[str, str],
        *,
        k: int = 15,
        max_kmer_owners: int | None = None,
    ) -> None:
        if k < 1:
            raise ValueError("k must be at least 1")
        if max_kmer_owners is not None and max_kmer_owners < 1:
            raise ValueError("max_kmer_owners must be positive")
        if not references:
            raise ValueError("at least one reference is required")

        normalized_by_id: dict[str, str] = {}
        for raw_id, raw_sequence in references.items():
            reference_id = str(raw_id).strip()
            if not reference_id:
                raise ValueError("reference identifiers must not be empty")
            if reference_id in normalized_by_id:
                raise ValueError(f"duplicate reference identifier: {reference_id!r}")
            normalized_by_id[reference_id] = normalize_sequence(raw_sequence)

        ids_by_sequence: dict[str, list[str]] = defaultdict(list)
        for reference_id, reference in normalized_by_id.items():
            ids_by_sequence[reference].append(reference_id)
        groups = tuple(
            sorted(
                ((tuple(sorted(ids)), sequence) for sequence, ids in ids_by_sequence.items()),
                key=lambda item: item[0],
            )
        )
        owners: dict[str, set[tuple[str, ...]]] = defaultdict(set)
        for aliases, reference in groups:
            for kmer in canonical_unique_kmers(reference, k):
                owners[kmer].add(aliases)

        self.k = k
        self.max_kmer_owners = max_kmer_owners
        self.references = tuple(groups)
        self._sequence_by_aliases = {aliases: sequence for aliases, sequence in groups}
        self._owners = {
            kmer: tuple(sorted(kmer_owners))
            for kmer, kmer_owners in owners.items()
            if max_kmer_owners is None or len(kmer_owners) <= max_kmer_owners
        }

    def shortlist(
        self,
        sequence: str,
        *,
        top_n: int = 5,
        min_kmer_score: float = 0.0,
    ) -> tuple[CandidateScore, ...]:
        """Return deterministic candidates using distinct specificity-weighted k-mers."""

        if top_n < 1:
            raise ValueError("top_n must be at least 1")
        query_kmers = canonical_unique_kmers(sequence, self.k)
        scores: dict[tuple[str, ...], float] = defaultdict(float)
        hits: dict[tuple[str, ...], int] = defaultdict(int)
        for kmer in query_kmers:
            owners = self._owners.get(kmer, ())
            if not owners:
                continue
            weight = 1.0 / len(owners)
            for aliases in owners:
                scores[aliases] += weight
                hits[aliases] += 1
        ranked = sorted(scores, key=lambda aliases: (-scores[aliases], aliases))
        return tuple(
            CandidateScore(aliases, scores[aliases], hits[aliases])
            for aliases in ranked[:top_n]
            if scores[aliases] >= min_kmer_score
        )

    def sequence_for(self, aliases: tuple[str, ...]) -> str:
        """Return the normalized sequence for an alias group."""

        return self._sequence_by_aliases[aliases]

    @property
    def alias_groups(self) -> tuple[tuple[str, ...], ...]:
        return tuple(aliases for aliases, _ in self.references)


def _align_candidate(
    query: str,
    reference: str,
    candidate: CandidateScore,
) -> AlignmentEvidence:
    # Align the shorter complete sequence inside the longer one. This exposes
    # truncation through coverage while avoiding an arbitrary terminal-gap cost.
    if len(query) <= len(reference):
        result = edlib.align(
            query,
            reference,
            mode="HW",
            task="locations",
            additionalEqualities=edlib_iupac_equalities(),
        )
        start, end = min(
            (int(start), int(end))
            for start, end in result["locations"]
            if start is not None and end is not None
        )
        query_start, query_end = 0, len(query)
        reference_start, reference_end = start, end + 1
    else:
        result = edlib.align(
            reference,
            query,
            mode="HW",
            task="locations",
            additionalEqualities=edlib_iupac_equalities(),
        )
        start, end = min(
            (int(start), int(end))
            for start, end in result["locations"]
            if start is not None and end is not None
        )
        query_start, query_end = start, end + 1
        reference_start, reference_end = 0, len(reference)
    edit_distance = int(result["editDistance"])
    query_aligned = query_end - query_start
    reference_aligned = reference_end - reference_start
    denominator = max(query_aligned, reference_aligned, 1)
    identity = max(0.0, 1.0 - edit_distance / denominator)
    return AlignmentEvidence(
        aliases=candidate.aliases,
        kmer_score=candidate.kmer_score,
        matching_kmers=candidate.matching_kmers,
        edit_distance=edit_distance,
        identity=identity,
        query_coverage=query_aligned / len(query),
        reference_coverage=reference_aligned / len(reference),
        query_start=query_start,
        query_end=query_end,
        reference_start=reference_start,
        reference_end=reference_end,
    )


def assign_sequence(
    sequence: str,
    index: ReferenceIndex,
    *,
    left_motif: str | None = None,
    right_motif: str | None = None,
    motif_max_edits: int = 2,
    min_insert_length: int = 0,
    max_insert_length: int | None = None,
    expected_insert_length: int | None = None,
    top_n: int = 5,
    min_kmer_score: float = 1.0,
    min_identity: float = 0.90,
    min_query_coverage: float = 0.90,
    min_reference_coverage: float = 0.90,
    min_identity_margin: float = 0.02,
    fallback_align_all: bool = True,
) -> AssignmentCall:
    """Assign one sequence using k-mer candidates followed by edlib validation.

    Supplying motifs makes successful extraction mandatory. Missing motifs are
    never silently replaced with full-read classification.
    """

    if (left_motif is None) != (right_motif is None):
        raise ValueError("left_motif and right_motif must be supplied together")
    for name, value in (
        ("min_identity", min_identity),
        ("min_query_coverage", min_query_coverage),
        ("min_reference_coverage", min_reference_coverage),
        ("min_identity_margin", min_identity_margin),
    ):
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must be between 0 and 1")
    read = normalize_sequence(sequence)
    extraction: InsertExtraction | None = None
    orientation: Literal["forward", "reverse", "unknown"] = "forward"
    query = read
    if left_motif is not None and right_motif is not None:
        extraction = extract_insert(
            read,
            left_motif,
            right_motif,
            max_edits=motif_max_edits,
            min_insert_length=min_insert_length,
            max_insert_length=max_insert_length,
            expected_insert_length=expected_insert_length,
        )
        if extraction.status == "motif_missing":
            return AssignmentCall(
                "motif_missing", (), "unknown", None, None, None, None,
                extraction.reason, (), (), extraction,
            )
        if extraction.status == "ambiguous":
            return AssignmentCall(
                "ambiguous", (), "unknown", None, None, None, None,
                extraction.reason, (), (), extraction,
            )
        assert extraction.sequence is not None
        query = extraction.sequence
        orientation = extraction.orientation
    if not query:
        return AssignmentCall(
            "no_match", (), orientation, query, None, None, None,
            "the extracted query is empty", (), (), extraction,
        )

    candidates = index.shortlist(query, top_n=top_n, min_kmer_score=min_kmer_score)
    if not candidates and fallback_align_all:
        candidates = tuple(CandidateScore(aliases, 0.0, 0) for aliases in index.alias_groups)
    if not candidates:
        return AssignmentCall(
            "no_match", (), orientation, query, None, None, None,
            "no reference met the k-mer evidence threshold", (), (), extraction,
        )

    alignments = tuple(
        sorted(
            (
                _align_candidate(query, index.sequence_for(candidate.aliases), candidate)
                for candidate in candidates
            ),
            key=lambda hit: (
                -hit.identity,
                -min(hit.query_coverage, hit.reference_coverage),
                hit.edit_distance,
                -hit.kmer_score,
                hit.aliases,
            ),
        )
    )
    best = alignments[0]
    second = alignments[1] if len(alignments) > 1 else None
    margin = best.identity - second.identity if second is not None else None
    passes = (
        best.identity >= min_identity
        and best.query_coverage >= min_query_coverage
        and best.reference_coverage >= min_reference_coverage
    )
    if not passes and fallback_align_all and len(candidates) < len(index.alias_groups):
        already_aligned = {candidate.aliases for candidate in candidates}
        fallback_candidates = tuple(
            CandidateScore(aliases, 0.0, 0)
            for aliases in index.alias_groups
            if aliases not in already_aligned
        )
        candidates = (*candidates, *fallback_candidates)
        alignments = tuple(
            sorted(
                (
                    *alignments,
                    *(
                        _align_candidate(
                            query,
                            index.sequence_for(candidate.aliases),
                            candidate,
                        )
                        for candidate in fallback_candidates
                    ),
                ),
                key=lambda hit: (
                    -hit.identity,
                    -min(hit.query_coverage, hit.reference_coverage),
                    hit.edit_distance,
                    -hit.kmer_score,
                    hit.aliases,
                ),
            )
        )
        best = alignments[0]
        second = alignments[1] if len(alignments) > 1 else None
        margin = best.identity - second.identity if second is not None else None
        passes = (
            best.identity >= min_identity
            and best.query_coverage >= min_query_coverage
            and best.reference_coverage >= min_reference_coverage
        )
    if not passes:
        return AssignmentCall(
            "no_match", (), orientation, query, best, second, margin,
            "best alignment did not meet identity and coverage thresholds",
            candidates, alignments, extraction,
        )
    second_is_plausible = second is not None and (
        second.identity >= min_identity
        and second.query_coverage >= min_query_coverage
        and second.reference_coverage >= min_reference_coverage
    )
    if second_is_plausible and margin is not None and margin < min_identity_margin:
        return AssignmentCall(
            "ambiguous", (), orientation, query, best, second, margin,
            "best and second alignments do not meet the identity margin",
            candidates, alignments, extraction,
        )
    status: AssignmentStatus = (
        "assigned_alias_set" if len(best.aliases) > 1 else "assigned_unique"
    )
    return AssignmentCall(
        status, best.aliases, orientation, query, best, second, margin,
        "best alignment met identity, coverage, and margin thresholds",
        candidates, alignments, extraction,
    )
