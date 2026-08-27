"""Full-length references: insert designs joined to their constant flanks.

An oligo-pool tool emits the variable insert of each design and nothing else,
but the molecule that gets sequenced is the whole amplicon: a constant 5' region
carrying the primer site, RBS and start codon, then the insert, then a constant
3' region carrying the fusion partner, tag, terminator and the other primer
site.  Consensus over the insert alone therefore reports nothing about the ~900
constant bases that were sequenced anyway.

This module turns insert references into full-length ones.  The constant regions
can be given directly, or derived from a single example construct - one binder
already inserted into the vector - by finding which insert it contains and
taking what lies either side.  Deriving them is the friendlier route: the
example is a file the bench already has, whereas two long strings pasted into a
config are two long strings to get wrong.

The region the pipeline reconstructs is bounded by, and excludes, the two primer
anchors, which is what :func:`nanopore3.assignment.extract_insert` returns for a
motif pair.  A full-length reference is therefore the *inner* flanks around the
insert, and :meth:`Flanks.flank` is the one place that convention is written
down.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

import edlib

from .sequence import normalize_sequence, reverse_complement

if TYPE_CHECKING:  # pragma: no cover - typing only
    from collections.abc import Iterable, Sequence

# Long enough to be unique in a plasmid-scale amplicon, short enough to survive
# nanopore error at the read ends where quality is worst.
DEFAULT_ANCHOR_LENGTH = 20

# A derived flank set is only trustworthy if the example construct really does
# contain one of the inserts.
DEFAULT_MINIMUM_TEMPLATE_IDENTITY = 0.90


class FlankError(ValueError):
    """A flank definition that would produce misleading references."""


@dataclass(frozen=True, slots=True)
class Flanks:
    """The constant regions either side of an insert, and the primer anchors.

    ``upstream`` and ``downstream`` are the full constant regions including the
    anchors; ``left_anchor`` and ``right_anchor`` are the primer sites that bound
    the reconstructed region and are *excluded* from it.
    """

    upstream: str
    downstream: str
    left_anchor: str
    right_anchor: str
    # True when the flanks were read *out of* the references, so those references
    # already carry them. Such a library must be trimmed to the primer anchors,
    # never flanked again - doing both produced references with the backbone
    # duplicated at each end.
    references_include_flanks: bool = False

    def transform(self, sequence: str) -> str:
        """Put one reference record into the coordinate system the pipeline uses.

        Every reference ends up as ``inner_upstream + insert + inner_downstream``:
        the region a read spans, bounded by the primer anchors and excluding them.
        A library supplied as inserts gets there by having the constant regions
        joined on; a library supplied already assembled gets there by having the
        anchors trimmed off.  Same destination, opposite operations, which is why
        the distinction has to be carried rather than inferred.
        """

        return self.trim(sequence) if self.references_include_flanks else self.flank(sequence)

    def trim(self, assembled: str) -> str:
        """Strip the primer anchors from an already-assembled reference."""

        sequence = normalize_sequence(assembled)
        start, end = len(self.left_anchor), len(sequence) - len(self.right_anchor)
        if end <= start:
            raise FlankError(
                f"reference of {len(sequence)} nt is shorter than its two "
                f"{len(self.left_anchor)} nt primer anchors"
            )
        return sequence[start:end]

    @property
    def inner_upstream(self) -> str:
        """The 5' constant sequence inside the anchor."""

        return self.upstream[len(self.left_anchor) :]

    @property
    def inner_downstream(self) -> str:
        """The 3' constant sequence inside the anchor."""

        return self.downstream[: len(self.downstream) - len(self.right_anchor)]

    def flank(self, insert: str) -> str:
        """Return the full-length reference for one insert."""

        return self.inner_upstream + normalize_sequence(insert) + self.inner_downstream

    def insert_span(self, full_length: int) -> tuple[int, int]:
        """Half-open coordinates of the insert inside a full-length reference."""

        start = len(self.inner_upstream)
        return start, full_length - len(self.inner_downstream)

    def insert_anchors(self) -> tuple[str, str]:
        """The motif pair that bounds the insert, for insert-scoped analyses.

        Chimera detection profiles a read window by window against the reference
        set, which only discriminates where the references differ.  In a
        full-length library most of every reference is identical, so that work is
        done on the insert alone: these are the inner ends of the constant
        regions, the same boundaries an insert-only run would use.
        """

        length = len(self.left_anchor)
        left = self.inner_upstream[-length:]
        right = self.inner_downstream[:length]
        if not left or not right:
            raise FlankError(
                "the constant regions are too short to yield insert anchors of "
                f"{length} nt inside the primer anchors"
            )
        both = self.upstream + self.downstream
        _require_unique_anchor(left, both, "insert 5'")
        _require_unique_anchor(right, both, "insert 3'")
        return left, right

    @property
    def constant_bases(self) -> int:
        """Constant bases added to every reference."""

        return len(self.inner_upstream) + len(self.inner_downstream)


def _require_unique_anchor(anchor: str, sequence: str, side: str) -> None:
    """An anchor that recurs could bound the wrong region."""

    if sequence.count(anchor) > 1:
        raise FlankError(
            f"the {side} anchor {anchor!r} occurs {sequence.count(anchor)} times in the "
            "constant regions; extraction could bound the wrong region. Use a longer "
            "anchor_length or trim the flank"
        )


def from_sequences(
    upstream: str,
    downstream: str,
    *,
    anchor_length: int = DEFAULT_ANCHOR_LENGTH,
) -> Flanks:
    """Build flanks from the two constant regions, anchors taken from their ends."""

    up = normalize_sequence(upstream)
    down = normalize_sequence(downstream)
    if anchor_length < 8:
        raise FlankError("anchor_length must be at least 8 to be specific")
    if len(up) < anchor_length or len(down) < anchor_length:
        raise FlankError(
            f"constant regions must each be at least anchor_length ({anchor_length}) nt; "
            f"got {len(up)} and {len(down)}"
        )
    left, right = up[:anchor_length], down[-anchor_length:]
    both = up + down
    _require_unique_anchor(left, both, "5'")
    _require_unique_anchor(right, both, "3'")
    return Flanks(upstream=up, downstream=down, left_anchor=left, right_anchor=right)


def from_template(
    template: str,
    inserts: Iterable[str],
    *,
    anchor_length: int = DEFAULT_ANCHOR_LENGTH,
    minimum_identity: float = DEFAULT_MINIMUM_TEMPLATE_IDENTITY,
) -> tuple[Flanks, str]:
    """Derive flanks from one example construct containing one of the inserts.

    Returns the flanks and the insert that matched, so the caller can report
    which design the example carried.  The template may be given in either
    orientation.
    """

    sequence = normalize_sequence(template)
    candidates = [normalize_sequence(insert) for insert in inserts]
    if not sequence or not candidates:
        raise FlankError("a template sequence and at least one insert are required")

    best: tuple[float, int, int, str] | None = None
    for oriented in (sequence, reverse_complement(sequence)):
        for insert in candidates:
            if not insert or len(insert) > len(oriented):
                continue
            result = edlib.align(insert, oriented, mode="HW", task="locations")
            distance = result["editDistance"]
            if distance < 0:
                continue
            identity = 1.0 - distance / len(insert)
            start, end = result["locations"][0]
            key = (identity, -start, len(insert), insert)
            if best is None or key > best[:4]:
                best = (identity, -start, len(insert), insert)
                found = (oriented, start, end, insert, identity)

    if best is None or found[4] < minimum_identity:
        seen = f"{found[4]:.3f}" if best is not None else "none"
        raise FlankError(
            "no insert aligns to the template above "
            f"{minimum_identity:.2f} identity (best {seen}); the template is probably "
            "not a construct from this library"
        )

    oriented, start, end, insert, _identity = found
    upstream, downstream = oriented[:start], oriented[end + 1 :]
    if not upstream or not downstream:
        raise FlankError(
            "the matched insert reaches an end of the template, so one constant region "
            "would be empty; supply a template that includes both primer sites"
        )
    return from_sequences(upstream, downstream, anchor_length=anchor_length), insert


def from_assembled(
    references: Sequence[str],
    *,
    anchor_length: int = DEFAULT_ANCHOR_LENGTH,
    minimum_constant: int = 60,
) -> Flanks:
    """Derive the constant regions from references that are already full length.

    The third way to describe a whole-vector library, alongside naming the flanks
    and giving a template: supply the assembled constructs and let the shared
    backbone declare itself.  Every reference begins with the same 5' constant
    region and ends with the same 3' one, so the longest common prefix and suffix
    across the set *are* the flanks, and what remains is the insert.

    This needs no alignment - the constructs share exact sequence at both ends - and
    it is what makes per-region accuracy, insert-scoped thresholds and
    insert-scoped chimera detection available to a library that was never
    expressed as inserts-plus-backbone.

    A set that does not share a backbone produces a tiny prefix and suffix by
    chance, so ``minimum_constant`` refuses that case rather than reporting an
    insert span that is mostly vector.
    """

    sequences = [normalize_sequence(reference) for reference in references]
    sequences = [sequence for sequence in sequences if sequence]
    if len(sequences) < 2:
        raise FlankError(
            "deriving flanks from assembled references needs at least two of them; "
            "with one there is nothing to hold constant. Give flanks.upstream and "
            "flanks.downstream, or flanks.template, instead"
        )

    shortest = min(len(sequence) for sequence in sequences)
    first = sequences[0]
    prefix = 0
    while prefix < shortest and all(s[prefix] == first[prefix] for s in sequences):
        prefix += 1
    suffix = 0
    while (
        suffix < shortest - prefix
        and all(s[len(s) - 1 - suffix] == first[len(first) - 1 - suffix] for s in sequences)
    ):
        suffix += 1

    # Order matters: a set of identical references has a huge prefix and a zero
    # suffix, which would otherwise be reported as a missing backbone.
    if prefix + suffix >= shortest:
        raise FlankError(
            "the references have no variable region between their shared ends: "
            "they are identical over the whole of the shortest one"
        )
    if prefix < minimum_constant or suffix < minimum_constant:
        raise FlankError(
            f"the references share only {prefix} nt at the 5' end and {suffix} nt at "
            f"the 3' end, below the {minimum_constant} nt expected of a common "
            "backbone. Either they are not all the same construct, or they are "
            "insert-only - in which case name the flanks explicitly instead"
        )
    derived = from_sequences(
        first[:prefix], first[len(first) - suffix :], anchor_length=anchor_length
    )
    return replace(derived, references_include_flanks=True)


def flank_sequences(
    inserts: Sequence[tuple[str, str]],
    flanks: Flanks,
) -> list[tuple[str, str]]:
    """Flank named inserts, checking that none already carries the constant region."""

    built: list[tuple[str, str]] = []
    for name, insert in inserts:
        normalized = normalize_sequence(insert)
        if flanks.inner_upstream and flanks.inner_upstream[-12:] in normalized:
            raise FlankError(
                f"insert {name!r} already contains the 5' constant region; flanking it "
                "again would double that sequence"
            )
        built.append((name, flanks.flank(normalized)))
    return built


__all__ = [
    "DEFAULT_ANCHOR_LENGTH",
    "FlankError",
    "Flanks",
    "flank_sequences",
    "from_sequences",
    "from_template",
]
