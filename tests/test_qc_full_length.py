"""QC for a consensus that spans the whole amplicon, not just the insert."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3 import qc

# 30 nt of 5' vector, a 30 nt insert, 30 nt of 3' vector.
FIVE = "GGGTTTAAACCCGGGTTTAAACCCGGGAAA"
INSERT = "ATGAAACGCGTGGAGGAAAAGGTCAAAGAG"
THREE = "TTTCCCGGGAAATTTCCCGGGAAATTTCCC"
REFERENCE = FIVE + INSERT + THREE
SPANS = {
    "flank_5p": (0, len(FIVE)),
    "insert": (len(FIVE), len(FIVE) + len(INSERT)),
    "flank_3p": (len(FIVE) + len(INSERT), len(REFERENCE)),
}


def edits_by_region(consensus: str) -> dict[str, int]:
    return {r.name: r.edit_distance for r in qc.region_metrics(consensus, REFERENCE, SPANS)}


class RegionMetricsTests(unittest.TestCase):
    """Each edit must be charged to the region it falls in."""

    def test_a_perfect_consensus_has_no_edits_anywhere(self) -> None:
        self.assertEqual(edits_by_region(REFERENCE), {"flank_5p": 0, "insert": 0, "flank_3p": 0})

    def test_a_substitution_is_charged_to_its_own_region(self) -> None:
        for region, position in (
            ("flank_5p", 10),
            ("insert", len(FIVE) + 10),
            ("flank_3p", len(FIVE) + len(INSERT) + 10),
        ):
            mutated = REFERENCE[:position] + ("C" if REFERENCE[position] != "C" else "G")
            mutated += REFERENCE[position + 1 :]
            with self.subTest(region=region):
                counts = edits_by_region(mutated)
                self.assertEqual(counts[region], 1)
                self.assertEqual(sum(counts.values()), 1)

    def test_a_deletion_is_charged_to_its_own_region(self) -> None:
        position = len(FIVE) + 10
        counts = edits_by_region(REFERENCE[:position] + REFERENCE[position + 1 :])
        self.assertEqual(counts["insert"], 1)
        self.assertEqual(sum(counts.values()), 1)

    def test_an_insertion_is_charged_to_the_region_it_lands_in(self) -> None:
        """Pins edlib's CIGAR convention rather than assuming it."""

        position = len(FIVE) + 10
        counts = edits_by_region(REFERENCE[:position] + "AA" + REFERENCE[position:])
        self.assertEqual(sum(counts.values()), 2)
        self.assertEqual(counts["insert"], 2)

    def test_identity_is_per_region_not_diluted_by_the_others(self) -> None:
        position = len(FIVE) + 10
        # Pick a base that differs from the original, or nothing is mutated.
        replacement = "C" if REFERENCE[position] != "C" else "G"
        mutated = REFERENCE[:position] + replacement + REFERENCE[position + 1 :]
        regions = {r.name: r for r in qc.region_metrics(mutated, REFERENCE, SPANS)}
        self.assertAlmostEqual(regions["insert"].identity, 1 - 1 / len(INSERT))
        self.assertEqual(regions["flank_5p"].identity, 1.0)

    def test_no_spans_means_no_region_reporting(self) -> None:
        self.assertEqual(qc.region_metrics(REFERENCE, REFERENCE, {}), ())


class FullLengthCodingTests(unittest.TestCase):
    """The ORF is found inside the consensus instead of assembled around it."""

    START = "ATGCAGCTT"
    STOP = "CATCACCACCACCATCACTAA"

    def build(self, insert: str) -> str:
        return "GGGTTTAAACCCGGG" + self.START + insert + self.STOP + "TTTCCCGGGAAA"

    def test_an_in_frame_orf_with_no_internal_stop_passes(self) -> None:
        frame, stops, _protein, index, _detail = qc.coding_checks_in_consensus(
            self.build("AAAGAGGAAGCG"), start_anchor=self.START, stop_anchor=self.STOP
        )
        self.assertEqual((frame, stops), ("pass", "pass"))
        self.assertIsNone(index)

    def test_a_length_off_a_codon_boundary_is_a_frameshift(self) -> None:
        frame, stops, _protein, _index, detail = qc.coding_checks_in_consensus(
            self.build("AAAGAGGAAGC"), start_anchor=self.START, stop_anchor=self.STOP
        )
        self.assertEqual(frame, "fail")
        self.assertEqual(stops, "not_evaluable")
        self.assertIn("codon boundary", detail)

    def test_an_internal_stop_is_found_and_located(self) -> None:
        frame, stops, _protein, index, detail = qc.coding_checks_in_consensus(
            self.build("AAATAAGAAGCG"), start_anchor=self.START, stop_anchor=self.STOP
        )
        self.assertEqual((frame, stops), ("pass", "fail"))
        self.assertEqual(index, 4)
        self.assertIn("residue 5", detail)

    def test_a_missing_anchor_is_not_evaluable_rather_than_a_failure(self) -> None:
        frame, stops, _protein, _index, detail = qc.coding_checks_in_consensus(
            "ACGT" * 40, start_anchor=self.START, stop_anchor=self.STOP
        )
        self.assertEqual((frame, stops), ("not_evaluable", "not_evaluable"))
        self.assertIn("anchor not found", detail)

    def test_anchors_are_found_through_sequencing_error(self) -> None:
        consensus = self.build("AAAGAGGAAGCG").replace(self.START, "ATGCAGCTA", 1)
        frame, stops, _protein, _index, _detail = qc.coding_checks_in_consensus(
            consensus, start_anchor=self.START, stop_anchor=self.STOP
        )
        self.assertEqual((frame, stops), ("pass", "pass"))


class EvaluateConsensusModeTests(unittest.TestCase):
    """Both consensus shapes are supported, and they do not interfere."""

    START, STOP = "ATGCAGCTT", "CATCACCACCACCATCACTAA"

    def test_full_length_mode_evaluates_the_orf_inside_the_consensus(self) -> None:
        consensus = "GGGTTTAAACCCGGG" + self.START + "AAAGAGGAAGCG" + self.STOP + "TTTCCC"
        result = qc.evaluate_consensus(
            consensus, consensus, coding_anchors=(self.START, self.STOP)
        )
        self.assertEqual(result.overall, "pass")
        self.assertEqual(result.reading_frame, "pass")
        self.assertEqual(result.internal_stops, "pass")

    def test_insert_mode_still_assembles_the_orf_from_constants(self) -> None:
        result = qc.evaluate_consensus(
            "AAAGAGGAAGCG",
            "AAAGAGGAAGCG",
            upstream_constant=self.START,
            downstream_constant=self.STOP,
        )
        self.assertEqual(result.reading_frame, "pass")
        self.assertEqual(result.internal_stops, "pass")

    def test_without_either_the_coding_checks_stay_not_evaluable(self) -> None:
        result = qc.evaluate_consensus("ACGTACGT", "ACGTACGT")
        self.assertEqual(result.reading_frame, "not_evaluable")
        self.assertEqual(result.internal_stops, "not_evaluable")

    def test_region_columns_appear_in_the_row_only_when_spans_are_given(self) -> None:
        plain = qc.evaluate_consensus(REFERENCE, REFERENCE).to_dict()
        self.assertNotIn("insert_identity", plain)
        with_spans = qc.evaluate_consensus(REFERENCE, REFERENCE, spans=SPANS).to_dict()
        self.assertEqual(with_spans["insert_edit_distance"], 0)
        self.assertEqual(with_spans["insert_length"], len(INSERT))
        self.assertIn("flank_3p_identity", with_spans)


if __name__ == "__main__":
    unittest.main()
