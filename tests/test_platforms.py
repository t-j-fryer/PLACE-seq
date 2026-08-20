"""Cross-platform comparison: nanopore clone picking against pooled Illumina."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3 import platforms


def clone(
    block: int,
    key: str,
    grade: str,
    *,
    well="A01",
    plate="P1",
    enc="A",
    kind="consensus",
    seq=None,
):
    return platforms.Clone(
        culture_plate=plate,
        well_id=well,
        encoding=enc,
        block=block,
        design_key=key,
        grade=grade,
        kind=kind,
        sequence_id=seq if seq is not None else f"{key}-seq",
    )


def read(block: int, key: str, category: str, *, group="SUMO_A", enc="A", edit=0, length=300):
    return platforms.IlluminaRead(
        group=group,
        encoding=enc,
        block=block,
        expected_block=block,
        design_key=key,
        category=category,
        edit_distance=edit,
        reference_length=length,
    )


class ReferenceNameTests(unittest.TestCase):
    def test_both_platforms_separators_parse_to_the_same_design(self) -> None:
        """Illumina writes A|Block_9_x where we write A_Block_9_x."""

        self.assertEqual(
            platforms.parse_reference_name("A|Block_9_dTF090_x"),
            platforms.parse_reference_name("A_Block_9_dTF090_x"),
        )

    def test_an_unencoded_name_parses_with_an_empty_encoding(self) -> None:
        self.assertEqual(platforms.parse_reference_name("Block_15_dTF083_1"), ("", 15, "dTF083_1"))

    def test_an_unparseable_name_is_none_rather_than_a_guess(self) -> None:
        self.assertIsNone(platforms.parse_reference_name("chim-abc123"))


class DistinctSequencesPerWellTests(unittest.TestCase):
    def test_wells_that_yielded_nothing_are_counted(self) -> None:
        """Empty wells are the point of the figure, so they set the denominator."""

        histogram = platforms.distinct_sequences_per_well(
            [clone(1, "d1", "perfect")], culture_plates=1, wells_per_plate=95
        )
        self.assertEqual(histogram[1], 1)
        self.assertEqual(histogram[0], 94)
        self.assertEqual(sum(histogram.values()), 95)

    def test_two_clones_with_the_same_sequence_count_once(self) -> None:
        """"Distinct sequences" means content, not name."""

        clones = [
            clone(1, "d1", "perfect", seq="same"),
            clone(1, "d2", "perfect", seq="same"),
        ]
        histogram = platforms.distinct_sequences_per_well(
            clones, culture_plates=1, wells_per_plate=95
        )
        self.assertEqual(histogram[1], 1)

    def test_counts_above_the_maximum_fold_into_the_last_bucket(self) -> None:
        clones = [clone(1, f"d{i}", "perfect", seq=f"s{i}") for i in range(5)]
        histogram = platforms.distinct_sequences_per_well(
            clones, culture_plates=1, wells_per_plate=95, maximum=3
        )
        self.assertEqual(histogram[3], 1)


class RecoveryTests(unittest.TestCase):
    def test_the_best_observation_of_a_design_wins(self) -> None:
        clones = [
            clone(1, "d1", "frameshift", well="A01"),
            clone(1, "d1", "perfect", well="B01"),
            clone(1, "d1", "screenable", well="C01"),
        ]
        self.assertEqual(platforms.nanopore_recovery(clones), {"d1": platforms.PERFECT})

    def test_a_chimera_is_not_an_observation_of_its_parent(self) -> None:
        clones = [clone(1, "d1", "frameshift", kind="chimera")]
        self.assertEqual(platforms.nanopore_recovery(clones), {})

    def test_illumina_membership_follows_the_reference_not_the_well(self) -> None:
        """A read from an A well that matched a B reference is evidence about B.

        Filtering on the well's group instead let references from outside a
        library count towards it and pushed recovery above 100%.
        """

        universe = {("B", 1, "d1")}
        reads = [read(1, "d1", "perfect_match", group="SUMO_A", enc="B")]
        self.assertEqual(platforms.illumina_recovery(reads, universe), {"d1": platforms.PERFECT})
        self.assertEqual(platforms.illumina_recovery(reads, {("A", 1, "d1")}), {})

    def test_recovery_fractions_are_over_the_designed_library(self) -> None:
        fractions = platforms.recovery_fractions({"d1": platforms.PERFECT}, 4)
        self.assertEqual(fractions[platforms.PERFECT], 25.0)
        self.assertEqual(fractions[platforms.SCREENABLE], 0.0)

    def test_categories_map_onto_the_same_three_outcomes(self) -> None:
        self.assertEqual(platforms.GRADE_CLASS["screenable"], platforms.SCREENABLE)
        self.assertEqual(platforms.CATEGORY_CLASS["nonsynonymous"], platforms.SCREENABLE)
        self.assertEqual(platforms.CATEGORY_CLASS["frameshift"], platforms.OTHER)


class MatchedDepthTests(unittest.TestCase):
    def test_effort_is_wells_that_yielded_a_clone_of_the_block(self) -> None:
        clones = [
            clone(1, "d1", "perfect", well="A01"),
            clone(1, "d2", "perfect", well="A01"),  # same well, one unit of effort
            clone(1, "d3", "perfect", well="B01"),
            clone(2, "d4", "perfect", well="C01"),
        ]
        self.assertEqual(dict(platforms.wells_sequenced_per_block(clones)), {1: 2, 2: 1})

    def test_subsampling_is_deterministic_for_a_seed(self) -> None:
        reads = [read(1, f"d{i}", "perfect_match") for i in range(50)]
        depths = {("SUMO_A", 1): 10}
        first = platforms.subsample_by_block(reads, depths, seed=7)
        second = platforms.subsample_by_block(reads, depths, seed=7)
        self.assertEqual(len(first), 10)
        self.assertEqual([r.design_key for r in first], [r.design_key for r in second])

    def test_a_block_with_too_few_reads_yields_what_it_has(self) -> None:
        reads = [read(1, "d1", "perfect_match")]
        sample = platforms.subsample_by_block(reads, {("SUMO_A", 1): 10}, seed=7)
        self.assertEqual(len(sample), 1)

    def test_reads_are_pooled_by_the_well_block_not_the_matched_reference(self) -> None:
        """The sampling unit is the well's block; where the read landed is a finding."""

        stray = platforms.IlluminaRead(
            group="SUMO_A",
            encoding="A",
            block=17,
            expected_block=1,
            design_key="d1",
            category="perfect_match",
            edit_distance=0,
            reference_length=300,
        )
        sample = platforms.subsample_by_block([stray], {("SUMO_A", 1): 1}, seed=7)
        self.assertEqual(len(sample), 1)


class IdentityTests(unittest.TestCase):
    def test_illumina_identity_is_one_minus_edits_over_reference_length(self) -> None:
        reads = [read(1, "d1", "nonsynonymous", edit=3, length=300)]
        self.assertAlmostEqual(
            platforms.illumina_read_identities(reads, {("A", 1, "d1")})[0], 0.99
        )

    def test_a_reference_without_a_length_is_skipped_not_divided_by_zero(self) -> None:
        reads = [read(1, "d1", "perfect_match", length=0)]
        self.assertEqual(platforms.illumina_read_identities(reads, {("A", 1, "d1")}), [])

    def test_summary_reports_mean_median_and_exact_rate(self) -> None:
        summary = platforms.summarise_identities([1.0, 1.0, 0.9, 0.8])
        self.assertEqual(summary["n"], 4)
        self.assertAlmostEqual(summary["mean"], 0.925)
        self.assertAlmostEqual(summary["median"], 0.95)
        self.assertAlmostEqual(summary["exact_fraction"], 0.5)

    def test_an_empty_set_does_not_divide_by_zero(self) -> None:
        self.assertEqual(platforms.summarise_identities([])["mean"], 0.0)


if __name__ == "__main__":
    unittest.main()
