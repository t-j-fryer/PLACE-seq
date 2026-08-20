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

from dataclasses import dataclass
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
