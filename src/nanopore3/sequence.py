"""Portable sequence validation and normalization primitives.

The functions in this module are intentionally small and free of file-system or
process state.  They are safe to call inside workers created with the Windows
``spawn`` multiprocessing start method.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from functools import lru_cache

from .errors import Nanopore3Error

IUPAC_DNA_BASES = frozenset("ACGTRYSWKMBDHVN")

# Anything outside the unambiguous alphabet separates canonical k-mer runs.
_AMBIGUOUS_RUN = re.compile("[^ACGT]+")

_IUPAC_MEMBERS: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

_COMPLEMENT = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
    }
)


class SequenceValidationError(Nanopore3Error, ValueError):
    """Raised when a sequence contains unsupported symbols."""


def normalize_sequence(
    sequence: str,
    *,
    allow_empty: bool = False,
    convert_rna: bool = True,
) -> str:
    """Return an uppercase, whitespace-free IUPAC DNA sequence.

    Embedded spaces and line breaks are ignored so that wrapped FASTA content
    can be normalized directly.  RNA ``U`` is converted to ``T`` by default.
    Any other symbol is rejected rather than silently deleted.
    """

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be str, got {type(sequence).__name__}")
    normalized = "".join(sequence.split()).upper()
    if convert_rna:
        normalized = normalized.replace("U", "T")
    if not normalized and not allow_empty:
        raise SequenceValidationError("sequence is empty")
    invalid = sorted(set(normalized).difference(IUPAC_DNA_BASES))
    if invalid:
        rendered = ", ".join(repr(base) for base in invalid)
        raise SequenceValidationError(f"unsupported sequence symbol(s): {rendered}")
    return normalized


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of an IUPAC DNA/RNA sequence."""

    return normalize_sequence(sequence, allow_empty=True).translate(_COMPLEMENT)[::-1]


def iupac_compatible(left: str, right: str) -> bool:
    """Return whether two individual IUPAC symbols share a concrete base."""

    left_base = normalize_sequence(left)
    right_base = normalize_sequence(right)
    if len(left_base) != 1 or len(right_base) != 1:
        raise ValueError("iupac_compatible expects two single-base strings")
    return bool(_IUPAC_MEMBERS[left_base] & _IUPAC_MEMBERS[right_base])


@lru_cache(maxsize=1)
def edlib_iupac_equalities() -> tuple[tuple[str, str], ...]:
    """Return deterministic additional equalities for IUPAC-aware edlib calls."""

    symbols = sorted(IUPAC_DNA_BASES)
    return tuple(
        (left, right)
        for left in symbols
        for right in symbols
        if left != right and _IUPAC_MEMBERS[left] & _IUPAC_MEMBERS[right]
    )


def canonical_kmer(kmer: str) -> str:
    """Return the lexicographically canonical strand of an unambiguous k-mer."""

    normalized = normalize_sequence(kmer)
    if set(normalized).difference("ACGT"):
        raise SequenceValidationError("canonical k-mers must contain only A, C, G, and T")
    reverse = reverse_complement(normalized)
    return min(normalized, reverse)


def _iter_run_canonical_kmers(run: str, k: int) -> Iterator[str]:
    """Yield canonical k-mers of one unambiguous A/C/G/T run.

    The reverse complement of the whole run is built once, so each window costs
    two slices and one comparison rather than a per-window normalization and
    reverse complement.  ``run[start:start + k]`` reverse-complements to
    ``reverse[n - start - k:n - start]``, so the strand chosen here is identical
    to :func:`canonical_kmer` while the per-window validation is removed.
    """

    length = len(run)
    reverse = run.translate(_COMPLEMENT)[::-1]
    for start in range(length - k + 1):
        forward = run[start : start + k]
        complement = reverse[length - start - k : length - start]
        yield forward if forward < complement else complement


def canonical_unique_kmers(sequence: str, k: int) -> frozenset[str]:
    """Return distinct canonical A/C/G/T k-mers, skipping ambiguous windows."""

    if k < 1:
        raise ValueError("k must be at least 1")
    normalized = normalize_sequence(sequence, allow_empty=True)
    if len(normalized) < k:
        return frozenset()
    # An ambiguity code can never occur inside a canonical k-mer, so it simply
    # splits the sequence into independent unambiguous runs.
    kmers: set[str] = set()
    for run in _AMBIGUOUS_RUN.split(normalized):
        if len(run) >= k:
            kmers.update(_iter_run_canonical_kmers(run, k))
    return frozenset(kmers)

