"""Golden Gate fragment excision, assembly, and interchangeability."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from nanopore3.fragments import (
    FragmentError,
    FragmentLibrary,
    assemble_fragments,
    excised_candidates,
)

# GGTCTC + 1 spacer + 4 nt overhang + body + 4 nt overhang + 1 spacer + GAGACC,
# wrapped in primer padding that is discarded during digestion.
PAD_5, PAD_3 = "ACATAAGCGATCCCAAGGTC", "AGCTATAAGAATTGCCGGGC"


def oligo(overhang_5: str, body: str, overhang_3: str, pad: bool = True) -> str:
    core = f"GGTCTCA{overhang_5}{body}{overhang_3}TGAGACC"
    return f"{PAD_5}{core}{PAD_3}" if pad else core


class ExcisionTests(unittest.TestCase):
    def test_padding_and_recognition_sites_are_removed(self) -> None:
        candidates = excised_candidates(oligo("GCTT", "AAACCCGGGTTT", "AGTG"))
        self.assertIn("GCTT" + "AAACCCGGGTTT" + "AGTG", candidates)

    def test_an_oligo_without_a_bsai_pair_yields_nothing(self) -> None:
        self.assertEqual(excised_candidates("ACGTACGTACGTACGT"), ())

    def test_a_spurious_site_in_padding_produces_several_candidates(self) -> None:
        # A GGTCTC inside the primer region must not be assumed to be the real one.
        decoy = "GGTCTCAAAAA" + oligo("GCTT", "AAACCCGGGTTT", "AGTG", pad=False)
        self.assertGreater(len(excised_candidates(decoy)), 1)


class AssemblyTests(unittest.TestCase):
    def test_a_single_fragment_gene_drops_both_vector_overhangs(self) -> None:
        body = "AAACCCGGGTTTAAACCCGGG"
        gene = body
        fragments = assemble_fragments([oligo("GCTT", body, "AGTG")], gene)
        self.assertEqual(len(fragments), 1)
        self.assertEqual(fragments[0].overhang_pair, ("GCTT", "AGTG"))
        # The span starts before the gene because the 5' overhang is the vector's.
        self.assertEqual(fragments[0].start, -4)
        self.assertEqual(fragments[0].end, len(gene) + 4)

    def test_two_fragments_share_one_junction_overhang(self) -> None:
        left, right = "AAACCCGGGTTT", "TTTGGGCCCAAA"
        junction = "CATG"
        gene = left + junction + right
        fragments = assemble_fragments(
            [oligo("GCTT", left, junction), oligo(junction, right, "AGTG")], gene
        )
        self.assertEqual(len(fragments), 2)
        self.assertEqual(fragments[0].overhang_3, junction)
        self.assertEqual(fragments[1].overhang_5, junction)
        # The shared overhang appears once in the product, not twice.
        self.assertEqual(gene.count(junction), 1)
        self.assertEqual(fragments[1].start, len(left))

    def test_a_spurious_padding_site_does_not_break_reconstruction(self) -> None:
        body = "AAACCCGGGTTTAAACCCGGG"
        decoy = "GGTCTCAAAAA" + oligo("GCTT", body, "AGTG", pad=False)
        fragments = assemble_fragments([decoy], body)
        self.assertEqual(fragments[0].sequence, "GCTT" + body + "AGTG")

    def test_mismatched_junction_overhangs_are_rejected(self) -> None:
        left, right = "AAACCCGGGTTT", "TTTGGGCCCAAA"
        with self.assertRaises(FragmentError):
            assemble_fragments(
                [oligo("GCTT", left, "CATG"), oligo("GGGG", right, "AGTG")],
                left + "CATG" + right,
            )

    def test_a_gene_that_does_not_reconstruct_is_an_error(self) -> None:
        with self.assertRaises(FragmentError):
            assemble_fragments([oligo("GCTT", "AAACCCGGGTTT", "AGTG")], "TTTTTTTTTTTT")


_CSV = """Block,Sequence Name,Overhang1,Overhang2,DNA Fragment 1,DNA Fragment 2,Full Sequence
1,geneA,,,{a1},,{ga}
1,geneB,CATG,,{b1},{b2},{gb}
1,geneC,CATG,,{c1},{c2},{gc}
"""


class FragmentLibraryTests(unittest.TestCase):
    def setUp(self) -> None:
        left, altleft = "AAACCCGGGTTT", "AAACCCGGGAAA"
        right, altright = "TTTGGGCCCAAA", "TTTGGGCCCTTT"
        junction = "CATG"
        self.ga = "GGGCCCAAATTTGGGCCC"
        self.gb = left + junction + right
        self.gc = altleft + junction + altright
        text = _CSV.format(
            a1=oligo("GCTT", self.ga, "AGTG"),
            ga=self.ga,
            b1=oligo("GCTT", left, junction),
            b2=oligo(junction, right, "AGTG"),
            gb=self.gb,
            c1=oligo("GCTT", altleft, junction),
            c2=oligo(junction, altright, "AGTG"),
            gc=self.gc,
        )
        self._dir = tempfile.TemporaryDirectory()
        self.path = Path(self._dir.name) / "design.csv"
        self.path.write_text(text, encoding="utf-8")
        self.library = FragmentLibrary.from_full_info_csv(self.path)

    def tearDown(self) -> None:
        self._dir.cleanup()

    def test_every_gene_reconstructs_from_its_oligos(self) -> None:
        self.assertEqual(self.library.summary()["genes"], 3)
        self.assertEqual(self.library.summary()["fragments"], 5)

    def test_genes_are_looked_up_by_sequence_not_name(self) -> None:
        """Design tables and reference FASTAs use different identifiers."""

        self.assertEqual(self.library.gene_for_sequence(self.gb).name, "geneB")
        self.assertIsNone(self.library.gene_for_sequence("ACGTACGTACGT"))

    def test_fragments_sharing_an_overhang_pair_are_interchangeable(self) -> None:
        gene_b = self.library.gene_for_sequence(self.gb)
        # geneB and geneC share both junctions, so each slot has two variants:
        # exactly the mis-assembly a chimera is drawn from.
        for fragment in gene_b.fragments:
            variants = self.library.interchangeable(fragment)
            owners = {name for names in variants.values() for name in names}
            self.assertIn("geneB", owners)
            self.assertIn("geneC", owners)
        self.assertEqual(self.library.summary()["interchangeable_pairs"], 2)

    def test_a_unique_gene_shares_no_slot(self) -> None:
        gene_a = self.library.gene_for_sequence(self.ga)
        variants = self.library.interchangeable(gene_a.fragments[0])
        owners = {name for names in variants.values() for name in names}
        self.assertEqual(owners, {"geneA"})

    def test_a_missing_column_is_reported_clearly(self) -> None:
        bad = self.path.with_name("bad.csv")
        bad.write_text("Sequence Name,Full Sequence\nx,ACGT\n", encoding="utf-8")
        with self.assertRaisesRegex(FragmentError, "DNA Fragment"):
            FragmentLibrary.from_full_info_csv(bad)


if __name__ == "__main__":
    unittest.main()
