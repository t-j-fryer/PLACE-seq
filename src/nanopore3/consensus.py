"""Deterministic, portable reference-guided consensus construction.

The baseline implementation deliberately uses only edlib.  Optional POA backends can
be added behind the same result contract without making the portable installation
depend on platform-specific executables.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import math
import re
from typing import Iterable, Sequence

import edlib


_CIGAR_TOKEN = re.compile(r"(\d+)([=XID])")


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
    backend: str = "edlib_pileup"
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


def build_reference_consensus(
    reference: str,
    reads: Sequence[ConsensusRead],
    *,
    group_id: str,
    min_depth: int = 3,
    max_reads: int = 100,
    min_support: float = 0.60,
    seed: int = 0,
) -> ConsensusResult:
    """Build a deterministic consensus and retain complete contributor provenance."""

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
