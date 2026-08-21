"""The statistical base caller: an allele must beat the group's own error rate."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3.consensus import (
    ConsensusRead,
    binomial_upper_tail,
    build_reference_consensus,
)

REFERENCE = (
    "ATGAAACGCGTGGAGGAAAAGGTCAAAGAGATTTTCGAAAAGGAACTGGAGAAATTAATCAAAGAAGGC"
    "GTAAAAGAGGAAGAAGCGAAGAAAATCGTTGAAGAACGTAAAGAGAAGAACATTGAGTACTGG"
)


def reads(sequences: list[str], *, quality: int = 30) -> list[ConsensusRead]:
    return [
        ConsensusRead(f"r{i}", s, tuple([quality] * len(s)))
        for i, s in enumerate(sequences)
    ]


def substitute(sequence: str, position: int, base: str) -> str:
    return sequence[:position] + base + sequence[position + 1 :]


class BinomialTailTests(unittest.TestCase):
    def test_matches_hand_computed_values(self) -> None:
        # P(X >= 1) for n=1, p=0.5
        self.assertAlmostEqual(binomial_upper_tail(1, 1, 0.5), 0.5)
        # P(X >= 2) for n=2, p=0.5 is 0.25
        self.assertAlmostEqual(binomial_upper_tail(2, 2, 0.5), 0.25)
        # P(X >= 1) for n=2, p=0.5 is 0.75
        self.assertAlmostEqual(binomial_upper_tail(1, 2, 0.5), 0.75)

    def test_edges_do_not_explode(self) -> None:
        self.assertEqual(binomial_upper_tail(0, 10, 0.1), 1.0)
        self.assertEqual(binomial_upper_tail(11, 10, 0.1), 0.0)

    def test_a_large_count_at_a_small_rate_is_vanishing(self) -> None:
        self.assertLess(binomial_upper_tail(30, 100, 0.02), 1e-12)


class StatisticalCallerTests(unittest.TestCase):
    def test_a_clean_clone_reproduces_its_reference(self) -> None:
        result = build_reference_consensus(
            REFERENCE, reads([REFERENCE] * 20), group_id="g", min_depth=6
        )
        self.assertEqual(result.sequence, REFERENCE)
        self.assertEqual(result.status, "consensus_pass")
        self.assertEqual(result.ambiguous_bases, 0)
        self.assertEqual(result.mixed_positions, 0)

    def test_a_real_substitution_in_every_read_is_called(self) -> None:
        mutant = substitute(REFERENCE, 40, "T" if REFERENCE[40] != "T" else "A")
        result = build_reference_consensus(
            REFERENCE, reads([mutant] * 20), group_id="g", min_depth=6
        )
        self.assertEqual(result.sequence, mutant)
        self.assertEqual(result.mixed_positions, 0)

    def test_a_two_allele_position_becomes_ambiguous_not_confident(self) -> None:
        """The case a flat 0.60 support fraction called a confident base."""

        other = "T" if REFERENCE[40] != "T" else "A"
        mutant = substitute(REFERENCE, 40, other)
        population = [mutant] * 19 + [REFERENCE] * 10  # 66% / 34%
        result = build_reference_consensus(
            REFERENCE, reads(population), group_id="g", min_depth=6
        )
        self.assertEqual(result.status, "mixed_variants")
        self.assertEqual(result.mixed_positions, 1)
        self.assertEqual(result.sequence[40], "N")

    def test_that_same_split_would_have_passed_the_old_support_rule(self) -> None:
        """Guards the reason for the change: 66% clears a 0.60 fraction."""

        other = "T" if REFERENCE[40] != "T" else "A"
        population = [substitute(REFERENCE, 40, other)] * 19 + [REFERENCE] * 10
        majority = build_reference_consensus(
            REFERENCE, reads(population), group_id="g", min_depth=6, caller="majority"
        )
        self.assertEqual(majority.sequence[40], other)
        self.assertEqual(majority.mixed_positions, 0)

    def test_one_stray_read_does_not_create_a_variant(self) -> None:
        other = "T" if REFERENCE[40] != "T" else "A"
        population = [REFERENCE] * 29 + [substitute(REFERENCE, 40, other)]
        result = build_reference_consensus(
            REFERENCE, reads(population), group_id="g", min_depth=6
        )
        self.assertEqual(result.sequence, REFERENCE)
        self.assertEqual(result.mixed_positions, 0)

    def test_low_quality_bases_are_not_counted(self) -> None:
        other = "T" if REFERENCE[40] != "T" else "A"
        good = [ConsensusRead(f"g{i}", REFERENCE, tuple([30] * len(REFERENCE))) for i in range(10)]
        junk = [
            ConsensusRead(
                f"b{i}", substitute(REFERENCE, 40, other), tuple([2] * len(REFERENCE))
            )
            for i in range(10)
        ]
        result = build_reference_consensus(
            REFERENCE, good + junk, group_id="g", min_depth=6, min_base_quality=10
        )
        self.assertEqual(result.sequence, REFERENCE)
        self.assertGreater(result.low_quality_bases, 0)

    def test_the_background_error_rate_is_reported(self) -> None:
        result = build_reference_consensus(
            REFERENCE, reads([REFERENCE] * 20), group_id="g", min_depth=6
        )
        self.assertGreater(result.background_error_rate, 0.0)
        self.assertLessEqual(result.background_error_rate, 0.25)

    def test_depth_below_the_floor_is_not_a_consensus(self) -> None:
        result = build_reference_consensus(
            REFERENCE, reads([REFERENCE] * 3), group_id="g", min_depth=6
        )
        self.assertEqual(result.status, "low_depth")

    def test_the_result_is_deterministic(self) -> None:
        population = [REFERENCE] * 15 + [substitute(REFERENCE, 12, "T")] * 15
        first = build_reference_consensus(
            REFERENCE, reads(population), group_id="g", min_depth=6, max_reads=20
        )
        second = build_reference_consensus(
            REFERENCE, reads(population), group_id="g", min_depth=6, max_reads=20
        )
        self.assertEqual(first.sequence, second.sequence)
        self.assertEqual(first.contributor_ids, second.contributor_ids)


if __name__ == "__main__":
    unittest.main()


class DeletionThresholdTests(unittest.TestCase):
    """Deletions are a majority decision at the configured support, not at 50%.

    Taking a deletion above 50% instead silently shortened 76 clones in a full
    run: nanopore deletion rates at awkward positions sit in the 50-60% band, so
    the gap between the two thresholds is exactly where the platform's error mode
    lives.
    """

    def population(self, deleted_fraction: float, total: int = 100) -> list[ConsensusRead]:
        deleted = round(total * deleted_fraction)
        with_gap = REFERENCE[:40] + REFERENCE[41:]
        return reads([with_gap] * deleted + [REFERENCE] * (total - deleted))

    def test_a_deletion_in_most_reads_but_under_the_floor_is_not_taken(self) -> None:
        result = build_reference_consensus(
            REFERENCE, self.population(0.55), group_id="g", min_depth=6, min_support=0.60
        )
        self.assertEqual(len(result.sequence), len(REFERENCE))

    def test_a_deletion_above_the_floor_is_taken(self) -> None:
        result = build_reference_consensus(
            REFERENCE, self.population(0.80), group_id="g", min_depth=6, min_support=0.60
        )
        self.assertEqual(len(result.sequence), len(REFERENCE) - 1)

    def test_a_deletion_minority_never_reports_mixture(self) -> None:
        """Deletion rates run at 20-25% per position with no mixture present."""

        result = build_reference_consensus(
            REFERENCE, self.population(0.25), group_id="g", min_depth=6
        )
        self.assertEqual(result.mixed_positions, 0)
        self.assertEqual(result.sequence, REFERENCE)

    def test_the_floor_follows_the_configured_support(self) -> None:
        loose = build_reference_consensus(
            REFERENCE, self.population(0.55), group_id="g", min_depth=6, min_support=0.50
        )
        self.assertEqual(len(loose.sequence), len(REFERENCE) - 1)


class DeletionDenominatorTests(unittest.TestCase):
    """The deletion share is measured against every read that spanned the position.

    Quality decides whether a base is counted; it must not decide how many reads
    were there. Excluding low-quality bases from the denominator raised the
    deletion share at positions sitting in the 55-65% band and shortened 71 clones
    in a full run.
    """

    def mixture(self, deleted: int, low_quality_bases: int, total: int = 20):
        with_gap = REFERENCE[:40] + REFERENCE[41:]
        out = [
            ConsensusRead(f"d{i}", with_gap, tuple([30] * len(with_gap)))
            for i in range(deleted)
        ]
        intact = total - deleted
        for i in range(intact):
            # Some reads carry the base at low quality: counted as spanning the
            # position, not counted as a vote.
            quality = 2 if i < low_quality_bases else 30
            qualities = [30] * len(REFERENCE)
            qualities[40] = quality
            out.append(ConsensusRead(f"k{i}", REFERENCE, tuple(qualities)))
        return out

    def test_low_quality_bases_still_count_as_reads_spanning(self) -> None:
        # 11 of 20 reads show the deletion: 55%, below the 0.60 floor, so the base
        # stays - and it must stay even when most of the other reads are low
        # quality at that position.
        result = build_reference_consensus(
            REFERENCE,
            self.mixture(deleted=11, low_quality_bases=8),
            group_id="g",
            min_depth=6,
            min_support=0.60,
        )
        self.assertEqual(len(result.sequence), len(REFERENCE))

    def test_a_dominant_deletion_is_still_taken(self) -> None:
        result = build_reference_consensus(
            REFERENCE,
            self.mixture(deleted=16, low_quality_bases=0),
            group_id="g",
            min_depth=6,
            min_support=0.60,
        )
        self.assertEqual(len(result.sequence), len(REFERENCE) - 1)
