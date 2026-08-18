"""Reference-free clustering of the reads in one well.

Two clones of the same construct can be 95% identical, differing only across a
short variable insert.  Nanopore read error (5-10%) is therefore *larger* than
the true difference between them, so clustering on overall read similarity
cannot work: the noise dominates the signal.

This module clusters on **co-varying positions** instead.  Reads are projected
onto a seed read's coordinates, positions where the reads systematically disagree
are identified, and each read is reduced to its alleles at those positions alone.
The invariant backbone contributes nothing either way, which is what makes the
method insensitive to how much sequence the clones share.

Random sequencing error is scattered independently across reads, so it produces
no consistent minor allele at any one position and rarely survives the frequency
floor.  A systematic error hotspot (a homopolymer, say) can survive it, so a
split is additionally required to be supported by several positions at once:
one position is an artefact, ten co-varying positions are a haplotype.

**Measured limit.**  Separation depends on the true divergence between clones
relative to the read error rate.  Clones sharing 2.4 kb of backbone and differing
across a 900 nt insert separate cleanly at 6% read error.  Clones differing across
only a 300 nt insert in 2.7 kb are ~94% identical -- a 6% divergence, equal to the
error rate -- and there the within-clone and between-clone distance distributions
touch and clusters can mix.  Treat a cluster as evidence, not proof, when the
variable region is a small fraction of the amplicon, and prefer the
reference-guided path where designed references exist.

``min_allele_fraction`` sets the smallest detectable minor clone: a clone present
at a lower fraction than this leaves no minor allele above the floor and is
absorbed into the majority, so it is a detection limit, not a bug.

The scaffold is a *polished* seed rather than a raw read.  A raw read carries its
own 5-10% errors, and every one of them looks like a variant position to every
other read, scattering spurious variants across the invariant backbone and
inflating within-clone distances until real clones stop separating.  One round of
majority polishing removes them.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

import edlib

from .consensus import ConsensusRead, build_reference_consensus

_CIGAR = re.compile(r"(\d+)([=XID])")
GAP = "-"
_REFINEMENT_ROUNDS = 5


@dataclass(frozen=True, slots=True)
class ReadCluster:
    """One putative clone: the reads supporting it and its consensus."""

    cluster_id: int
    read_ids: tuple[str, ...]
    consensus: str
    size: int
    variant_positions: tuple[int, ...]
    mean_mismatches_to_centre: float


@dataclass(frozen=True, slots=True)
class ClusteringResult:
    clusters: tuple[ReadCluster, ...]
    unassigned_read_ids: tuple[str, ...]
    variant_positions: tuple[int, ...]
    seed_read_id: str
    reason: str


def _project_all(reads, seed: str) -> dict[str, list[str]]:
    return {read_id: _project(sequence, seed) for read_id, sequence in reads}


def _majority(projected: dict[str, list[str]], read_ids: Sequence[str]) -> str:
    """Plain majority base per column, used only to polish the scaffold."""

    width = len(next(iter(projected.values())))
    bases = []
    for index in range(width):
        column = Counter(
            projected[read_id][index]
            for read_id in read_ids
            if projected[read_id][index] != GAP
        )
        if column:
            bases.append(max(sorted(column.items()), key=lambda item: item[1])[0])
    return "".join(bases)


def _project(query: str, seed: str) -> list[str]:
    """Return the query's base at each seed position, or ``GAP`` where absent."""

    result = edlib.align(query, seed, mode="NW", task="path")
    cigar = result.get("cigar") or ""
    bases = [GAP] * len(seed)
    qi = si = 0
    for count_text, operation in _CIGAR.findall(cigar):
        count = int(count_text)
        if operation in "=X":
            for _ in range(count):
                if si < len(seed):
                    bases[si] = query[qi]
                qi += 1
                si += 1
        elif operation == "I":       # extra base in the read
            qi += count
        elif operation == "D":       # base absent from the read
            si += count
    return bases


def choose_seed(reads: Sequence[tuple[str, str]]) -> int:
    """Index of the read nearest the median length.

    A median-length read is a safer scaffold than the longest, which is the one
    most likely to carry a chimeric or adapter-laden artefact.
    """

    lengths = sorted(len(sequence) for _, sequence in reads)
    target = lengths[len(lengths) // 2]
    return min(range(len(reads)), key=lambda i: (abs(len(reads[i][1]) - target), i))


def find_variant_positions(
    columns: Sequence[Counter[str]],
    *,
    min_depth: int,
    min_allele_fraction: float,
) -> list[int]:
    """Positions carrying a consistent minor allele, ignoring gaps."""

    positions = []
    for index, column in enumerate(columns):
        observed = {base: n for base, n in column.items() if base != GAP}
        depth = sum(observed.values())
        if depth < min_depth or len(observed) < 2:
            continue
        ranked = sorted(observed.values(), reverse=True)
        if ranked[1] / depth >= min_allele_fraction:
            positions.append(index)
    return positions


def _centre(vectors: dict[str, list[str]], read_ids: Sequence[str]) -> list[str]:
    """Per-position majority allele of a group, ignoring uncalled positions."""

    width = len(next(iter(vectors.values())))
    centre = []
    for index in range(width):
        column = Counter(
            vectors[read_id][index]
            for read_id in read_ids
            if vectors[read_id][index] != GAP
        )
        centre.append(max(sorted(column.items()), key=lambda item: item[1])[0] if column else GAP)
    return centre


def _distance(left: Sequence[str], right: Sequence[str]) -> tuple[int, int]:
    """Mismatches and comparable positions between two allele vectors."""

    mismatches = comparable = 0
    for a, b in zip(left, right, strict=True):
        if a == GAP or b == GAP:
            continue
        comparable += 1
        if a != b:
            mismatches += 1
    return mismatches, comparable


def cluster_reads(
    reads: Sequence[tuple[str, str]],
    *,
    min_cluster_size: int = 6,
    min_allele_fraction: float = 0.20,
    max_mismatch_fraction: float = 0.40,
    min_variant_positions: int = 3,
    min_comparable_fraction: float = 0.5,
    min_support: float = 0.60,
) -> ClusteringResult:
    """Group reads from one well into putative clones without a reference.

    ``min_variant_positions`` is the guard against splitting on a single
    systematic error hotspot: fewer co-varying positions than this and the reads
    are treated as one clone.

    ``min_comparable_fraction`` is the guard against a noisy distance.  An indel
    near a variant position leaves that position uncalled, and two reads that
    overlap at only a handful of positions produce a mismatch ratio driven by
    chance.  Without this floor the within-clone and between-clone distance
    distributions overlap and no threshold separates them.
    """

    if min_cluster_size < 1:
        raise ValueError("min_cluster_size must be positive")
    if not 0.0 < min_allele_fraction < 0.5:
        raise ValueError("min_allele_fraction must be between 0 and 0.5")
    if len(reads) < min_cluster_size:
        return ClusteringResult((), tuple(rid for rid, _ in reads), (), "", "too few reads")

    seed_index = choose_seed(reads)
    seed_id, raw_seed = reads[seed_index]
    # Polish the scaffold before looking for variants: a raw seed's own errors
    # would otherwise appear as variant positions to every other read.
    draft = _project_all(reads, raw_seed)
    seed = _majority(draft, [read_id for read_id, _ in reads])
    projected = _project_all(reads, seed)

    columns: list[Counter[str]] = [Counter() for _ in seed]
    for bases in projected.values():
        for index, base in enumerate(bases):
            columns[index][base] += 1
    positions = find_variant_positions(
        columns,
        min_depth=max(min_cluster_size, 2),
        min_allele_fraction=min_allele_fraction,
    )

    if len(positions) < min_variant_positions:
        # Not enough structured disagreement to justify a split: one clone.
        consensus = _cluster_consensus(
            reads, [read_id for read_id, _ in reads], seed, min_support
        )
        return ClusteringResult(
            (
                ReadCluster(
                    0, tuple(read_id for read_id, _ in reads), consensus,
                    len(reads), tuple(positions), 0.0,
                ),
            ),
            (), tuple(positions), seed_id,
            f"{len(positions)} variant position(s), fewer than the {min_variant_positions} "
            "required to split",
        )

    vectors = {
        read_id: [bases[index] for index in positions]
        for read_id, bases in projected.items()
    }
    # Seed clusters from the most common complete haplotypes, then absorb the
    # rest by nearest centre. Ordering by frequency keeps the result stable.
    order = sorted(
        vectors,
        key=lambda read_id: (
            -sum(1 for base in vectors[read_id] if base != GAP),
            read_id,
        ),
    )
    required = max(min_variant_positions, int(min_comparable_fraction * len(positions)))
    centres: list[list[str]] = []
    members: list[list[str]] = []
    for read_id in order:
        vector = vectors[read_id]
        best = None
        for index, centre in enumerate(centres):
            mismatches, comparable = _distance(vector, centre)
            if comparable < required:
                continue
            if mismatches / comparable <= max_mismatch_fraction:
                if best is None or mismatches < best[1]:
                    best = (index, mismatches)
        if best is None:
            centres.append(list(vector))
            members.append([read_id])
        else:
            members[best[0]].append(read_id)

    # A single greedy pass fixes each read against whichever centre existed when
    # it was seen, so an early centre can absorb a read from another clone.
    # Recomputing centres from their members and reassigning repairs that; a few
    # rounds are enough and it is run to a fixed point or a small cap.
    for _ in range(_REFINEMENT_ROUNDS):
        centres = [_centre(vectors, group) for group in members if group]
        regrouped: list[list[str]] = [[] for _ in centres]
        for read_id in order:
            vector = vectors[read_id]
            best = None
            for index, centre in enumerate(centres):
                mismatches, comparable = _distance(vector, centre)
                if comparable < required:
                    continue
                if mismatches / comparable <= max_mismatch_fraction:
                    if best is None or mismatches < best[1]:
                        best = (index, mismatches)
            if best is None:
                centres.append(list(vector))
                regrouped.append([read_id])
            else:
                regrouped[best[0]].append(read_id)
        if regrouped == members:
            break
        members = regrouped

    big = [group for group in members if len(group) >= min_cluster_size]
    if len(big) < 2:
        # One cluster reaching the floor is not evidence of several clones; it is
        # one clone whose reads scattered. Splitting here would silently discard
        # the reads that fell short, which reads as lower depth rather than as a
        # failed split.
        everything = [read_id for read_id, _ in reads]
        return ClusteringResult(
            (
                ReadCluster(
                    0, tuple(sorted(everything)),
                    _cluster_consensus(reads, everything, seed, min_support),
                    len(everything), tuple(positions), 0.0,
                ),
            ),
            (), tuple(positions), seed_id,
            f"{len(positions)} variant positions but only {len(big)} cluster(s) "
            f"reached {min_cluster_size} reads: treated as one clone",
        )

    clusters: list[ReadCluster] = []
    unassigned: list[str] = []
    for index, group in sorted(enumerate(members), key=lambda item: -len(item[1])):
        if len(group) < min_cluster_size:
            unassigned.extend(group)
            continue
        distances = [_distance(vectors[r], centres[index]) for r in group]
        mean_mismatch = (
            sum(m for m, _ in distances) / len(distances) if distances else 0.0
        )
        clusters.append(
            ReadCluster(
                len(clusters),
                tuple(sorted(group)),
                _cluster_consensus(reads, group, seed, min_support),
                len(group),
                tuple(positions),
                mean_mismatch,
            )
        )
    return ClusteringResult(
        tuple(clusters), tuple(sorted(unassigned)), tuple(positions), seed_id,
        f"{len(clusters)} cluster(s) over {len(positions)} variant positions",
    )


def _cluster_consensus(
    reads: Sequence[tuple[str, str]],
    read_ids: Sequence[str],
    scaffold: str,
    min_support: float,
) -> str:
    """Consensus for one cluster, scaffolded on the polished seed.

    This delegates to the same pileup the reference-guided path uses, so
    insertions are handled and the result is not capped by the projection used
    for clustering, which deliberately ignores them.
    """

    wanted = set(read_ids)
    selected = [
        ConsensusRead(read_id, sequence)
        for read_id, sequence in reads
        if read_id in wanted
    ]
    if not selected:
        return ""
    result = build_reference_consensus(
        scaffold, selected, group_id="cluster",
        min_depth=1, max_reads=len(selected), min_support=min_support,
    )
    return result.sequence


def _consensus_from(
    projected: dict[str, list[str]],
    read_ids: Sequence[str],
    columns: Sequence[Counter[str]] | None,
    min_support: float,
) -> str:
    """Majority-vote consensus over the seed's coordinates for one group."""

    if columns is None:
        width = len(next(iter(projected.values())))
        columns = [Counter() for _ in range(width)]
        for read_id in read_ids:
            for index, base in enumerate(projected[read_id]):
                columns[index][base] += 1
    bases = []
    for column in columns:
        observed = {base: n for base, n in column.items() if base != GAP}
        depth = sum(observed.values())
        if not depth:
            continue
        base, count = max(sorted(observed.items()), key=lambda item: item[1])
        # Below the support floor the position is genuinely unresolved, so it is
        # reported as N rather than silently taking a plurality winner.
        bases.append(base if count / depth >= min_support else "N")
    return "".join(bases)


__all__ = [
    "ClusteringResult",
    "ReadCluster",
    "choose_seed",
    "cluster_reads",
    "find_variant_positions",
]
