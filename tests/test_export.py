"""Graded, browsable consensus output organised by plate and well."""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from nanopore3.export import (
    GRADES,
    grade_consensus,
    normalized_well,
    safe_name,
    write_consensus_tree,
)


def consensus(**overrides):
    row = {
        "consensus_id": "cons-1", "sample_id": "S", "plate_id": "RP06",
        "well_id": "A1", "reference_library_id": "lib", "reference_ids": "Block_1_gene",
        "status": "consensus_pass", "culture_plate": "CP_A", "n_reads_available": "30",
        "n_reads_used": "30", "mean_depth": "30.0", "min_depth": "28",
        "ambiguous_bases": "0", "backend": "portable", "sequence_sha256": "",
        "failure_reason": "",
    }
    row.update(overrides)
    return row


def qc(**overrides):
    row = {
        "consensus_id": "cons-1", "full_amplicon": "pass", "expected_length": "pass",
        "reading_frame": "pass", "internal_stops": "pass", "overall": "pass",
        "alignment_edit_distance": "0", "alignment_identity": "1.0",
        "reason": "all evaluable criteria passed", "protein_length": "100",
    }
    row.update(overrides)
    return row


class GradeTests(unittest.TestCase):
    def test_an_exact_match_is_perfect(self) -> None:
        self.assertEqual(grade_consensus(consensus(), qc()), "perfect")

    def test_substitutions_but_passing_qc_are_screenable(self) -> None:
        self.assertEqual(
            grade_consensus(consensus(), qc(alignment_edit_distance="4")), "screenable"
        )

    def test_low_depth_wins_over_everything(self) -> None:
        self.assertEqual(
            grade_consensus(
                consensus(status="low_depth"), qc(overall="not_evaluable")
            ),
            "low_depth",
        )

    def test_a_mixed_well_is_reported_as_heterogeneous(self) -> None:
        """Ambiguity codes explain the identity failure that follows from them."""

        self.assertEqual(
            grade_consensus(
                consensus(ambiguous_bases="7"), qc(full_amplicon="fail", overall="fail")
            ),
            "heterogeneous",
        )

    def test_a_broken_reading_frame_outranks_length_and_identity(self) -> None:
        self.assertEqual(
            grade_consensus(
                consensus(),
                qc(reading_frame="fail", expected_length="fail",
                   full_amplicon="fail", internal_stops="not_evaluable", overall="fail"),
            ),
            "frameshift",
        )

    def test_a_premature_stop_outranks_length(self) -> None:
        self.assertEqual(
            grade_consensus(
                consensus(),
                qc(internal_stops="fail", expected_length="fail", overall="fail"),
            ),
            "premature_stop",
        )

    def test_every_grade_produced_is_a_declared_grade(self) -> None:
        cases = [
            (consensus(), qc()),
            (consensus(status="low_depth"), qc(overall="not_evaluable")),
            (consensus(ambiguous_bases="1"), qc()),
            (consensus(), qc(expected_length="fail", overall="fail")),
            (consensus(), qc(full_amplicon="fail", overall="fail")),
        ]
        for row, qc_row in cases:
            self.assertIn(grade_consensus(row, qc_row), GRADES)


class NamingTests(unittest.TestCase):
    def test_wells_are_padded_so_a_listing_sorts_correctly(self) -> None:
        self.assertEqual(normalized_well("A1"), "A01")
        self.assertEqual(normalized_well("h12"), "H12")
        self.assertEqual(normalized_well("A01"), "A01")

    def test_unsafe_characters_are_replaced(self) -> None:
        self.assertEqual(safe_name("a/b c:d"), "a_b_c_d")

    def test_long_names_shorten_but_stay_distinct(self) -> None:
        left, right = safe_name("x" * 200 + "A"), safe_name("x" * 200 + "B")
        self.assertLessEqual(len(left), 56)
        self.assertNotEqual(left, right)


class TreeTests(unittest.TestCase):
    def setUp(self) -> None:
        self._dir = tempfile.TemporaryDirectory()
        self.root = Path(self._dir.name) / "tree"

    def tearDown(self) -> None:
        self._dir.cleanup()

    def test_files_land_under_plate_and_well_with_grade_in_the_name(self) -> None:
        rows = [
            consensus(consensus_id="c1", well_id="A1"),
            consensus(consensus_id="c2", well_id="A10", reference_ids="Block_2_other"),
        ]
        qcs = {
            "c1": qc(consensus_id="c1"),
            "c2": qc(consensus_id="c2", alignment_edit_distance="3"),
        }
        summary = write_consensus_tree(rows, qcs, {"c1": "ACGT" * 20, "c2": "TTTT" * 20}, self.root)
        self.assertEqual(summary["files_written"], 2)
        first = self.root / "RP06" / "A01" / "RP06_A01__Block_1_gene__perfect.fasta"
        second = self.root / "RP06" / "A10" / "RP06_A10__Block_2_other__screenable.fasta"
        self.assertTrue(first.exists())
        self.assertTrue(second.exists())
        self.assertIn("grade=perfect", first.read_text())
        self.assertIn("culture_plate=CP_A", first.read_text())

    def test_a_consensus_without_a_sequence_is_indexed_but_writes_no_file(self) -> None:
        """A missing file must never have to be read as an oversight."""

        rows = [consensus(consensus_id="c1", status="low_depth", n_reads_used="2")]
        summary = write_consensus_tree(
            rows, {"c1": qc(consensus_id="c1", overall="not_evaluable")}, {}, self.root
        )
        self.assertEqual(summary["files_written"], 0)
        with (self.root / "index.csv").open() as handle:
            index = list(csv.DictReader(handle))
        self.assertEqual(len(index), 1)
        self.assertEqual(index[0]["grade"], "low_depth")
        self.assertEqual(index[0]["file"], "")

    def test_summary_reports_grades_overall_and_per_plate(self) -> None:
        rows = [
            consensus(consensus_id="c1", plate_id="RP06"),
            consensus(consensus_id="c2", plate_id="RP07", well_id="B2"),
        ]
        qcs = {"c1": qc(consensus_id="c1"), "c2": qc(consensus_id="c2")}
        summary = write_consensus_tree(rows, qcs, {"c1": "ACGT", "c2": "ACGT"}, self.root)
        self.assertEqual(summary["grades"]["perfect"], 2)
        self.assertEqual(summary["by_plate"]["RP06"]["perfect"], 1)
        self.assertEqual(summary["by_plate"]["RP07"]["perfect"], 1)
        written = json.loads((self.root / "summary.json").read_text())
        self.assertIn("perfect", written["descriptions"])

    def test_two_designs_in_one_well_both_get_files(self) -> None:
        """Polyclonal wells are the normal case here, not an error."""

        rows = [
            consensus(consensus_id="c1", reference_ids="Block_1_a"),
            consensus(consensus_id="c2", reference_ids="Block_2_b"),
        ]
        qcs = {"c1": qc(consensus_id="c1"), "c2": qc(consensus_id="c2")}
        write_consensus_tree(rows, qcs, {"c1": "AC", "c2": "GT"}, self.root)
        self.assertEqual(len(list((self.root / "RP06" / "A01").iterdir())), 2)


if __name__ == "__main__":
    unittest.main()
