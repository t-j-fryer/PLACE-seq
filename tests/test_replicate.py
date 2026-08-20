"""Concordance between a dedicated plate and the same plate inside a pool."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3 import replicate


def call(identity: str, sequence: str, *, kind: str = "consensus") -> replicate.Call:
    return replicate.Call(
        identity=identity,
        kind=kind,
        status="consensus_pass",
        reads=30,
        sequence=sequence,
    )


class WellNameTests(unittest.TestCase):
    def test_zero_padding_does_not_split_a_well(self) -> None:
        """The consensus table writes A1 and the export tree writes A01."""

        self.assertEqual(replicate._normalise_well("A1"), replicate._normalise_well("A01"))
        self.assertEqual(replicate._normalise_well("h12"), "H12")


class CompareWellsTests(unittest.TestCase):
    def test_same_clone_same_sequence_is_exact(self) -> None:
        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGT")}},
            {"A01": {"d1": call("d1", "ACGT")}},
        )
        self.assertEqual([o.category for o in outcomes], [replicate.MATCHED_EXACT])
        self.assertEqual(outcomes[0].exact, ("d1",))

    def test_same_clone_different_sequence_is_reported_with_a_distance(self) -> None:
        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGTACGT")}},
            {"A01": {"d1": call("d1", "ACGTACGA")}},
        )
        self.assertEqual(outcomes[0].category, replicate.MATCHED_DIFFERENT_SEQUENCE)
        self.assertEqual(outcomes[0].edit_distances, (("d1", 1),))

    def test_a_polyclonal_well_is_exact_only_if_every_shared_clone_is(self) -> None:
        """One matching clone must not certify a well with a second, differing one."""

        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGT"), "d2": call("d2", "TTTT")}},
            {"A01": {"d1": call("d1", "ACGT"), "d2": call("d2", "TTTA")}},
        )
        self.assertEqual(outcomes[0].category, replicate.MATCHED_DIFFERENT_SEQUENCE)
        self.assertEqual(outcomes[0].exact, ("d1",))
        self.assertEqual(outcomes[0].differing, ("d2",))

    def test_disjoint_clone_names_are_a_disagreement_not_a_match(self) -> None:
        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGT")}},
            {"A01": {"d2": call("d2", "ACGT")}},
        )
        self.assertEqual(outcomes[0].category, replicate.DIFFERENT_DESIGN)

    def test_each_side_alone_is_labelled_by_which_side_saw_it(self) -> None:
        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGT")}},
            {"B02": {"d2": call("d2", "ACGT")}},
        )
        self.assertEqual(
            {o.well_id: o.category for o in outcomes},
            {"A01": replicate.DEDICATED_ONLY, "B02": replicate.POOL_ONLY},
        )

    def test_a_missing_sequence_is_never_called_exact(self) -> None:
        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "")}},
            {"A01": {"d1": call("d1", "")}},
        )
        self.assertEqual(outcomes[0].category, replicate.MATCHED_DIFFERENT_SEQUENCE)
        self.assertEqual(outcomes[0].edit_distances, ())


class SummariseTests(unittest.TestCase):
    def test_concordance_is_measured_over_wells_the_dedicated_run_saw(self) -> None:
        """A well only the pool recovered is extra yield, not a disagreement."""

        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGT")}},
            {
                "A01": {"d1": call("d1", "ACGT")},
                "B01": {"d2": call("d2", "ACGT")},
                "B02": {"d3": call("d3", "ACGT")},
            },
        )
        summary = replicate.summarise(outcomes)
        self.assertEqual(summary["wells_total"], 3)
        self.assertEqual(summary["wells_dedicated"], 1)
        self.assertEqual(summary["wells_matched"], 1)
        self.assertEqual(summary["matched_fraction"], 1.0)
        self.assertEqual(summary["counts"][replicate.POOL_ONLY], 2)

    def test_clone_level_counts_are_reported_alongside_well_level(self) -> None:
        outcomes = replicate.compare_wells(
            {"A01": {"d1": call("d1", "ACGT"), "d2": call("d2", "TTTT")}},
            {"A01": {"d1": call("d1", "ACGT"), "d2": call("d2", "TTTA")}},
        )
        summary = replicate.summarise(outcomes)
        self.assertEqual(summary["clones_shared"], 2)
        self.assertEqual(summary["clones_exact"], 1)
        self.assertAlmostEqual(summary["exact_fraction"], 0.5)

    def test_matched_chimeras_are_counted(self) -> None:
        """A chimera's parent signature identifies it as a design name does."""

        parents = "Block_1_dTF001 >> Block_1_dTF002"
        outcomes = replicate.compare_wells(
            {"A01": {parents: call(parents, "ACGT", kind="chimera")}},
            {"A01": {parents: call(parents, "ACGT", kind="chimera")}},
        )
        self.assertEqual(replicate.summarise(outcomes)["chimeras_matched"], 1)

    def test_empty_input_does_not_divide_by_zero(self) -> None:
        summary = replicate.summarise([])
        self.assertEqual(summary["matched_fraction"], 0.0)
        self.assertEqual(summary["exact_fraction"], 0.0)


class FastaTests(unittest.TestCase):
    def test_multiline_records_are_joined_and_keyed_by_first_token(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "x.fasta"
            path.write_text(">a plate=RP01\nACGT\nTTTT\n>b\nGG\n", encoding="utf-8")
            self.assertEqual(replicate.read_fasta(path), {"a": "ACGTTTTT", "b": "GG"})

    def test_a_missing_file_is_empty_not_an_error(self) -> None:
        self.assertEqual(replicate.read_fasta(Path("/nonexistent/x.fasta")), {})


if __name__ == "__main__":
    unittest.main()
