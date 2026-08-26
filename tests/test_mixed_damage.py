"""The damage-signature grade, and what it is allowed to claim.

A mixed clone is not one thing.  These tests pin the distinction the grade rests
on: every contested position carrying G:C -> T:A is the 8-oxoguanine signature and
grades as a usable clone; anything else stays an unexplained mixture.
"""

from __future__ import annotations

import unittest

from nanopore3 import platforms
from nanopore3.export import GRADES, grade_consensus
from nanopore3.pipeline import parse_mixed_alleles, reading_frame_span
from nanopore3.qc import (
    classify_mixed_positions,
    designed_allele_fraction,
    mixed_signature,
    worst_protein_effect,
)

# ATG, seven codons, TAA: MKPGFWHD*. The TGG codon spans positions 15-17.
REFERENCE = "ATG" + "AAACCCGGGTTTTGGCATGAC" + "TAA"


class MixedSignatureTests(unittest.TestCase):
    def test_g_to_t_and_c_to_a_are_the_only_damage_signature(self) -> None:
        oxidative = classify_mixed_positions(
            ((9, "G", "GT"), (4, "C", "AC")), REFERENCE
        )
        self.assertTrue(all(p.oxidative for p in oxidative))
        self.assertEqual(mixed_signature(oxidative), "oxidative")

        for position, base, alleles in (
            (9, "G", "AG"),  # transition, the nanopore background spectrum
            (9, "G", "CG"),  # a transversion, but not the oxidative one
            (4, "C", "CT"),
        ):
            other = classify_mixed_positions(((position, base, alleles),), REFERENCE)
            self.assertFalse(other[0].oxidative, f"{base}>{alleles}")
            self.assertEqual(mixed_signature(other), "mixed")

    def test_one_unexplained_position_disqualifies_the_whole_clone(self) -> None:
        mixed = classify_mixed_positions(((9, "G", "GT"), (12, "T", "CT")), REFERENCE)
        self.assertEqual(mixed_signature(mixed), "mixed")

    def test_no_mixed_positions_has_no_signature(self) -> None:
        self.assertEqual(mixed_signature(classify_mixed_positions((), REFERENCE)), "")


class LocationTests(unittest.TestCase):
    def test_region_and_protein_effect_are_reported(self) -> None:
        spans = {"insert": (0, 12), "flank_3p": (12, len(REFERENCE))}
        mixed = classify_mixed_positions(
            ((9, "G", "GT"), (18, "G", "GT")),
            REFERENCE,
            spans=spans,
            reading_frame=(0, len(REFERENCE)),
        )
        self.assertEqual(mixed[0].region, "insert")
        self.assertEqual(mixed[1].region, "flank_3p")
        self.assertTrue(all(p.protein_effect for p in mixed))

    def test_a_position_outside_the_reading_frame_has_no_protein_effect(self) -> None:
        # The reading frame stops at 12, so a mixture past it cannot change the
        # protein - the case that matters least and must be visibly different.
        mixed = classify_mixed_positions(
            ((18, "G", "GT"),), REFERENCE, reading_frame=(0, 12)
        )
        self.assertEqual(mixed[0].protein_effect, "")
        self.assertIn(":", mixed[0].describe())

    def test_a_silent_mixture_is_distinguished_from_a_nonsense_one(self) -> None:
        # The TGG codon spans 15-17, so G>A at 17 makes TGA, a stop.
        nonsense = classify_mixed_positions(
            ((17, "G", "AG"),), REFERENCE, reading_frame=(0, len(REFERENCE))
        )
        self.assertEqual(nonsense[0].protein_effect, "nonsense")
        # The third base of the proline codon at 6-8 is silent either way.
        silent = classify_mixed_positions(
            ((8, "C", "AC"),), REFERENCE, reading_frame=(0, len(REFERENCE))
        )
        self.assertEqual(silent[0].protein_effect, "silent")

    def test_the_column_round_trips(self) -> None:
        written = "12:G>GT:0.6100,0.3900|880:C>AC:0.5000,0.5000"
        self.assertEqual(
            parse_mixed_alleles(written),
            ((12, "G", "GT", (0.61, 0.39)), (880, "C", "AC", (0.5, 0.5))),
        )
        self.assertEqual(parse_mixed_alleles(""), ())
        self.assertEqual(parse_mixed_alleles("nonsense"), ())

    def test_a_run_written_before_the_fractions_existed_still_parses(self) -> None:
        self.assertEqual(parse_mixed_alleles("12:G>GT"), ((12, "G", "GT", ()),))

    def test_the_reading_frame_is_located_not_assumed(self) -> None:
        reference = "TTTT" + REFERENCE + "GGGG"
        span = reading_frame_span(reference, REFERENCE[:6], REFERENCE[-6:])
        self.assertEqual(span, (4, 4 + len(REFERENCE)))


class DesignedAlleleTests(unittest.TestCase):
    """Whether the designed sequence is still in the well, which screening asks."""

    def test_the_designed_share_is_reported_per_position(self) -> None:
        mixed = classify_mixed_positions(
            ((9, "G", "GT", (0.62, 0.38)),), REFERENCE,
            reading_frame=(0, len(REFERENCE)),
        )
        self.assertAlmostEqual(mixed[0].designed_fraction, 0.62)
        self.assertEqual(mixed[0].variant_alleles, "T")
        self.assertAlmostEqual(designed_allele_fraction(mixed), 0.62)

    def test_a_clone_is_only_as_recoverable_as_its_weakest_position(self) -> None:
        mixed = classify_mixed_positions(
            # Fractions pair with the alleles in the order they are listed, so
            # "CT" with (0.79, 0.21) is C at 79% and the designed T at 21%.
            ((9, "G", "GT", (0.62, 0.38)), (12, "T", "CT", (0.79, 0.21))),
            REFERENCE,
        )
        # 0.21 for T at position 12, not the 0.62 at position 9.
        self.assertAlmostEqual(designed_allele_fraction(mixed), 0.21)

    def test_the_designed_base_can_be_absent_entirely(self) -> None:
        mixed = classify_mixed_positions(((9, "G", "AT", (0.55, 0.45)),), REFERENCE)
        self.assertEqual(mixed[0].designed_fraction, 0.0)

    def test_no_mixed_positions_means_no_answer_rather_than_zero(self) -> None:
        self.assertIsNone(designed_allele_fraction(()))

    def test_the_worst_effect_across_positions_is_reported(self) -> None:
        mixed = classify_mixed_positions(
            ((8, "C", "AC"), (17, "G", "AG")), REFERENCE,
            reading_frame=(0, len(REFERENCE)),
        )
        self.assertEqual({p.protein_effect for p in mixed}, {"silent", "nonsense"})
        self.assertEqual(worst_protein_effect(mixed), "nonsense")


class GradeTests(unittest.TestCase):
    consensus = {"status": "mixed_variants", "ambiguous_bases": "1"}

    def test_a_damage_signature_grades_above_an_unexplained_mixture(self) -> None:
        damage = grade_consensus(
            self.consensus, {"overall": "pass", "mixed_signature": "oxidative"}
        )
        self.assertEqual(damage, "mixed_damage")
        self.assertLess(GRADES.index("mixed_damage"), GRADES.index("mixed_variants"))

    def test_an_unexplained_mixture_is_unchanged(self) -> None:
        self.assertEqual(
            grade_consensus(self.consensus, {"overall": "pass", "mixed_signature": "mixed"}),
            "mixed_variants",
        )
        # A run whose QC predates the column must not silently become "damage".
        self.assertEqual(
            grade_consensus(self.consensus, {"overall": "pass"}), "mixed_variants"
        )

    def test_low_depth_still_wins_over_the_damage_signature(self) -> None:
        self.assertEqual(
            grade_consensus(
                {"status": "low_depth", "ambiguous_bases": "1"},
                {"overall": "pass", "mixed_signature": "oxidative"},
            ),
            "low_depth",
        )

    def test_a_damage_clone_counts_as_screenable_not_other(self) -> None:
        self.assertEqual(platforms.GRADE_CLASS["mixed_damage"], platforms.SCREENABLE)
        self.assertEqual(platforms.GRADE_CLASS["mixed_variants"], platforms.OTHER)


if __name__ == "__main__":
    unittest.main()
