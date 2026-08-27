"""Golden Gate fragment models for oligo-pool assemblies.

An oligo in a pool is not the sequence that ends up in the gene.  Each oligo
carries a BsaI recognition site and primer padding at both ends:

``[padding] GGTCTC N <4 nt overhang> ... <4 nt overhang> N GAGACC [padding]``

BsaI cuts outside its recognition site, so the retained portion runs from the
start of the 5' overhang to the end of the 3' overhang.  Adjacent fragments are
joined *through* a shared overhang, which therefore appears once in the product,
and the outermost overhangs are contributed by the vector.

Genes are assembled in **blocks**: sub-pools amplified with their own primer
pair and ligated in separate reactions.  Junction overhangs are designed to be
unique *within* a block and are deliberately reused *between* blocks, which is
safe because two blocks never share a tube.  Any interchangeability question is
therefore only meaningful within one block; comparing across blocks invents
recombinants that cannot physically exist.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from .errors import Nanopore3Error
from .sequence import normalize_sequence

BSAI_FORWARD = "GGTCTC"
BSAI_REVERSE = "GAGACC"
# BsaI cuts one base past its recognition site, leaving a four-base overhang.
_BSAI_SPACER = 1
OVERHANG_LENGTH = 4


class FragmentError(Nanopore3Error, ValueError):
    """A fragment table could not be parsed or did not reconstruct its gene."""


@dataclass(frozen=True, slots=True)
class Fragment:
    """One oligo's contribution to a gene, including both shared overhangs."""

    slot: int
    sequence: str
    start: int
    end: int

    @property
    def overhang_5(self) -> str:
        return self.sequence[:OVERHANG_LENGTH]

    @property
    def overhang_3(self) -> str:
        return self.sequence[-OVERHANG_LENGTH:]

    @property
    def overhang_pair(self) -> tuple[str, str]:
        """Fragments sharing this pair are interchangeable during assembly."""

        return (self.overhang_5, self.overhang_3)


@dataclass(frozen=True, slots=True)
class Gene:
    """A designed gene and the fragments it was assembled from."""

    name: str
    sequence: str
    fragments: tuple[Fragment, ...]
    # The sub-pool this gene was amplified and assembled in. Mis-assembly is
    # confined to one block, so this bounds every recombination hypothesis.
    block: str = ""


def _occurrences(text: str, motif: str) -> list[int]:
    found: list[int] = []
    index = text.find(motif)
    while index >= 0:
        found.append(index)
        index = text.find(motif, index + 1)
    return found


def excised_candidates(oligo: str) -> tuple[str, ...]:
    """Return every plausible post-digestion contribution of one oligo.

    Primer padding can itself contain a BsaI site, so the correct pair is the
    one that reconstructs the gene rather than simply the first found.
    """

    sequence = normalize_sequence(oligo, allow_empty=True)
    candidates = []
    for start in _occurrences(sequence, BSAI_FORWARD):
        body_start = start + len(BSAI_FORWARD) + _BSAI_SPACER
        for end in _occurrences(sequence, BSAI_REVERSE):
            body_end = end - _BSAI_SPACER
            if body_end - body_start >= 2 * OVERHANG_LENGTH:
                candidates.append(sequence[body_start:body_end])
    return tuple(candidates)


def assemble_fragments(oligos: list[str], gene: str) -> tuple[Fragment, ...]:
    """Choose one contribution per oligo that reconstructs ``gene`` exactly.

    Junction overhangs must match between neighbours, and the assembled product
    must equal the gene once the two vector overhangs are removed.  Raises when
    no combination reconstructs the gene, because a silently wrong fragment
    model would misreport every downstream assembly call.
    """

    options = [excised_candidates(oligo) for oligo in oligos]
    for slot, choices in enumerate(options):
        if not choices:
            raise FragmentError(f"fragment {slot + 1} contains no usable BsaI pair")
    target = normalize_sequence(gene)
    solution: list[Fragment] | None = None

    def walk(slot: int, assembled: str, chosen: list[Fragment]) -> bool:
        nonlocal solution
        if slot == len(options):
            # The outermost overhangs come from the vector, not from the gene.
            if assembled[OVERHANG_LENGTH:-OVERHANG_LENGTH] == target:
                solution = chosen
                return True
            return False
        for candidate in options[slot]:
            if slot == 0:
                start = 0
                merged = candidate
            elif assembled[-OVERHANG_LENGTH:] == candidate[:OVERHANG_LENGTH]:
                start = len(assembled) - OVERHANG_LENGTH
                merged = assembled + candidate[OVERHANG_LENGTH:]
            else:
                continue
            fragment = Fragment(slot, candidate, start, start + len(candidate))
            if walk(slot + 1, merged, [*chosen, fragment]):
                return True
        return False

    walk(0, "", [])
    if solution is None:
        raise FragmentError("no BsaI selection reconstructs the supplied gene")
    # Re-express spans in gene coordinates by dropping the 5' vector overhang.
    return tuple(
        Fragment(
            fragment.slot,
            fragment.sequence,
            fragment.start - OVERHANG_LENGTH,
            fragment.end - OVERHANG_LENGTH,
        )
        for fragment in solution
    )


class FragmentLibrary:
    """Genes, their fragments, and the interchangeable variants per overhang pair."""

    def __init__(self, genes: tuple[Gene, ...]) -> None:
        if not genes:
            raise FragmentError("a fragment library requires at least one gene")
        self.genes = genes
        self._by_sequence: dict[str, Gene] = {}
        for gene in genes:
            self._by_sequence.setdefault(gene.sequence, gene)
        self._by_block: dict[str, list[Gene]] = defaultdict(list)
        for gene in genes:
            self._by_block[gene.block].append(gene)
        # Keyed by (block, overhang pair): a shared pair only permits
        # mis-assembly when both genes were in the same reaction.
        variants: dict[tuple[str, tuple[str, str]], dict[str, list[str]]] = defaultdict(
            lambda: defaultdict(list)
        )
        for gene in genes:
            for fragment in gene.fragments:
                key = (gene.block, fragment.overhang_pair)
                variants[key][fragment.sequence].append(gene.name)
        self.variants_by_overhang: dict[
            tuple[str, tuple[str, str]], dict[str, tuple[str, ...]]
        ] = {
            key: {
                sequence: tuple(sorted(owners))
                for sequence, owners in sorted(by_sequence.items())
            }
            for key, by_sequence in sorted(variants.items())
        }

    @classmethod
    def from_full_info_csv(cls, path: str | Path) -> FragmentLibrary:
        """Load an oPool ``*_FULL_INFO.csv`` design table."""

        source = Path(path).expanduser().resolve(strict=True)
        genes: list[Gene] = []
        with source.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"Sequence Name", "Full Sequence"}
            missing = required.difference(reader.fieldnames or ())
            if missing:
                raise FragmentError(
                    f"{source} is missing column(s): {', '.join(sorted(missing))}"
                )
            columns = [
                name
                for name in (reader.fieldnames or ())
                if name.startswith("DNA Fragment")
            ]
            if not columns:
                raise FragmentError(f"{source} has no 'DNA Fragment' columns")
            for number, row in enumerate(reader, start=2):
                name = (row.get("Sequence Name") or "").strip()
                gene = (row.get("Full Sequence") or "").strip()
                if not name or not gene:
                    continue
                oligos = [row[column].strip() for column in columns if row[column].strip()]
                try:
                    fragments = assemble_fragments(oligos, gene)
                except (FragmentError, ValueError) as exc:
                    raise FragmentError(f"{source}:{number} ({name}): {exc}") from exc
                genes.append(
                    Gene(
                        name,
                        normalize_sequence(gene),
                        fragments,
                        block=(row.get("Block") or "").strip(),
                    )
                )
        return cls(tuple(genes))

    def gene_for_sequence(self, sequence: str) -> Gene | None:
        """Look up a gene by exact sequence.

        Design tables and reference FASTAs are joined by sequence because their
        identifiers differ between the two files.
        """

        return self._by_sequence.get(normalize_sequence(sequence, allow_empty=True))

    def interchangeable(self, block: str, fragment: Fragment) -> dict[str, tuple[str, ...]]:
        """Fragment variants that could occupy this slot *in the same block*.

        Overhangs repeat across blocks by design, so this must be scoped to one
        block. Scoping it globally reports recombinants between genes that were
        never in the same tube.
        """

        return self.variants_by_overhang.get((block, fragment.overhang_pair), {})

    def genes_in_block(self, block: str) -> tuple[Gene, ...]:
        """Every gene assembled in one sub-pool: the full recombination space.

        Junction overhangs are unique within a block, so a chimera cannot arise
        from a legitimate overhang swap. Real ones come from mis-ligation of
        similar overhangs or from PCR template switching, and both are bounded
        by the block rather than by the overhang pair.
        """

        return tuple(self._by_block.get(block, ()))

    @property
    def blocks(self) -> tuple[str, ...]:
        return tuple(sorted(self._by_block))

    def summary(self) -> dict[str, int]:
        pairs = self.variants_by_overhang
        return {
            "genes": len(self.genes),
            "blocks": len(self._by_block),
            "fragments": sum(len(gene.fragments) for gene in self.genes),
            "within_block_overhang_pairs": len(pairs),
            "shared_within_block": sum(1 for v in pairs.values() if len(v) > 1),
        }


__all__ = [
    "BSAI_FORWARD",
    "BSAI_REVERSE",
    "Fragment",
    "FragmentError",
    "FragmentLibrary",
    "Gene",
    "assemble_fragments",
    "excised_candidates",
]
