"""Positional reference matching: detect chimeric clones and keep their reads.

A chimeric molecule is a real clone in a real well, not an artefact to discard.
Assignment cannot describe one, because it asks a whole-read question -- "which
single reference is this?" -- and a chimera has no honest answer to it: an even
split falls below the identity floor and is rejected, while a lopsided one is
confidently attributed to whichever parent dominates.

This module asks a positional question instead.  Each read's insert is walked in
overlapping windows, and every window is matched to the reference library on its
own.  A clean clone yields the same reference in every window; a chimera yields a
run of one parent followed by a run of another.  Collapsing the per-window calls
into a signature therefore names both parents and locates the junction.

Comparing read *segments to references* rather than reads to each other is what
makes this work where read clustering does not: two designs sharing 94% of their
sequence are hard to separate read-by-read at nanopore error rates, but the
windows that differ still match their own parent unambiguously.

Reads sharing a signature are the same clone, so grouping by signature and
requiring the usual depth floor yields consensus sequences for chimeras on the
same terms as for anything else.  Signatures below the floor are noise and are
left alone.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

from .assignment import ReferenceIndex

DEFAULT_WINDOW = 90
DEFAULT_STEP = 45


@dataclass(frozen=True, slots=True)
class ReadSignature:
    """One read's positional match profile."""

    read_id: str
    windows: tuple[str | None, ...]
    signature: tuple[str, ...]

    @property
    def is_chimeric(self) -> bool:
        return len(self.signature) > 1

    @property
    def junction_windows(self) -> tuple[int, ...]:
        """Indices where the matched reference changes."""

        changes = []
        previous = None
        for index, name in enumerate(self.windows):
            if name is None:
                continue
            if previous is not None and name != previous:
                changes.append(index)
            previous = name
        return tuple(changes)


@dataclass(frozen=True, slots=True)
class ChimeraGroup:
    """Reads sharing one multi-parent signature: a putative chimeric clone."""

    signature: tuple[str, ...]
    read_ids: tuple[str, ...]
    size: int
    junction_window: int | None

    @property
    def label(self) -> str:
        return "__".join(self.signature)


def positional_profile(
    insert: str,
    index: ReferenceIndex,
    *,
    window: int = DEFAULT_WINDOW,
    step: int = DEFAULT_STEP,
    min_kmer_score: float = 1.0,
) -> tuple[str | None, ...]:
    """Best-matching reference for each sliding window of one insert."""

    if window < 1 or step < 1:
        raise ValueError("window and step must be positive")
    calls: list[str | None] = []
    for start in range(0, max(1, len(insert) - window + 1), step):
        candidates = index.shortlist(
            insert[start : start + window], top_n=1, min_kmer_score=min_kmer_score
        )
        calls.append(candidates[0].aliases[0] if candidates else None)
    return tuple(calls)


def collapse_signature(windows: Iterable[str | None]) -> tuple[str, ...]:
    """Reduce per-window calls to the ordered run of distinct references."""

    out: list[str] = []
    for name in windows:
        if name is None:
            continue
        if not out or out[-1] != name:
            out.append(name)
    return tuple(out)


def signature_for(
    read_id: str,
    insert: str,
    index: ReferenceIndex,
    *,
    window: int = DEFAULT_WINDOW,
    step: int = DEFAULT_STEP,
) -> ReadSignature:
    calls = positional_profile(insert, index, window=window, step=step)
    return ReadSignature(read_id, calls, collapse_signature(calls))


def group_chimeras(
    signatures: Sequence[ReadSignature],
    *,
    minimum_depth: int = 6,
) -> tuple[tuple[ChimeraGroup, ...], Counter[str]]:
    """Group chimeric reads by signature, keeping only groups deep enough.

    A single read can acquire a spurious extra segment from one noisy window, so
    a signature carried by fewer than ``minimum_depth`` reads is treated as noise
    rather than as a clone.  Real chimeras recur across every read of the clone.
    """

    by_signature: dict[tuple[str, ...], list[ReadSignature]] = defaultdict(list)
    outcome: Counter[str] = Counter()
    for item in signatures:
        if not item.signature:
            outcome["unmatched"] += 1
            continue
        if not item.is_chimeric:
            outcome["uniform"] += 1
            continue
        outcome["chimeric_read"] += 1
        by_signature[item.signature].append(item)

    groups = []
    for signature, members in by_signature.items():
        if len(members) < minimum_depth:
            outcome["below_depth"] += len(members)
            continue
        junctions = Counter(
            j for member in members for j in member.junction_windows
        )
        groups.append(
            ChimeraGroup(
                signature=signature,
                read_ids=tuple(sorted(m.read_id for m in members)),
                size=len(members),
                junction_window=junctions.most_common(1)[0][0] if junctions else None,
            )
        )
        outcome["in_group"] += len(members)
    groups.sort(key=lambda g: (-g.size, g.signature))
    return tuple(groups), outcome


def synthesise_reference(
    signature: Sequence[str],
    references: Mapping[str, str],
    junction_window: int | None,
    *,
    window: int = DEFAULT_WINDOW,
    step: int = DEFAULT_STEP,
) -> str:
    """Build a scaffold for a two-parent chimera by splicing at the junction.

    Only the two-parent case is spliced.  With three or more segments the
    junctions are not independently located well enough to trust a scaffold, so
    the caller should fall back to a reference-free consensus.
    """

    if len(signature) != 2 or junction_window is None:
        return ""
    left, right = references.get(signature[0], ""), references.get(signature[1], "")
    if not left or not right:
        return ""
    cut = min(junction_window * step, len(left))
    tail = min(cut, len(right))
    return left[:cut] + right[tail:]


__all__ = [
    "DEFAULT_STEP",
    "DEFAULT_WINDOW",
    "ChimeraGroup",
    "ReadSignature",
    "collapse_signature",
    "group_chimeras",
    "positional_profile",
    "signature_for",
    "synthesise_reference",
]
