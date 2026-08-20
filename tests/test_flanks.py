"""Full-length references built from insert designs and their constant flanks."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3 import flanks
from nanopore3.assignment import extract_insert

# Enough distinct sequence either side to be anchored and unique.
UPSTREAM = "CATAATCCGCACGCATCTGGTCGATCCCGCGAAATTAATACGACTCACTATAGGGAGAGCAGCTATGCAGCTT"
DOWNSTREAM = (
    "AGTGGATCCGAAACACCAGGCACGTCGGAAAGCGCTACCCCTGAATCAGTGTTTACCTTGGAGGACTTTGTG"
    "CATCACCACCACCATCACTAATGACTCGAGTCTGGTAAAGAAACCGCTGCTGCGAAATTTGGATTGGCGAATGGGACGC"
)
INSERT = "ATGAAACGCGTGGAGGAAAAGGTCAAAGAGATTTTCGAAAAGGAACTGGAGAAATTAATCAAAGAAGGC"


class BuildTests(unittest.TestCase):
    def test_a_full_length_reference_is_the_inner_flanks_around_the_insert(self) -> None:
        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        built = f.flank(INSERT)
        self.assertTrue(built.startswith(UPSTREAM[20:]))
        self.assertTrue(built.endswith(DOWNSTREAM[:-20]))
        self.assertEqual(len(built), len(UPSTREAM) - 20 + len(INSERT) + len(DOWNSTREAM) - 20)

    def test_the_insert_span_locates_the_variable_part(self) -> None:
        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        built = f.flank(INSERT)
        start, end = f.insert_span(len(built))
        self.assertEqual(built[start:end], INSERT)

    def test_a_short_anchor_is_refused_rather_than_matched_everywhere(self) -> None:
        with self.assertRaises(flanks.FlankError):
            flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=4)

    def test_a_flank_shorter_than_the_anchor_is_refused(self) -> None:
        with self.assertRaises(flanks.FlankError):
            flanks.from_sequences("ACGTACGTAC", DOWNSTREAM, anchor_length=20)

    def test_a_repeated_anchor_is_refused_because_it_could_bound_the_wrong_region(self) -> None:
        repeated = "AAACCCGGGTTTAAACCCGGGTTT"
        with self.assertRaises(flanks.FlankError) as raised:
            flanks.from_sequences(repeated + "ACGT", repeated + DOWNSTREAM, anchor_length=12)
        self.assertIn("occurs", str(raised.exception))

    def test_flanking_an_insert_that_already_carries_the_flank_is_refused(self) -> None:
        """Doubling the constant region would be silently wrong, not loud."""

        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        with self.assertRaises(flanks.FlankError):
            flanks.flank_sequences([("d1", f.flank(INSERT))], f)

    def test_named_inserts_keep_their_names(self) -> None:
        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        built = flanks.flank_sequences([("d1", INSERT), ("d2", INSERT[:30])], f)
        self.assertEqual([name for name, _ in built], ["d1", "d2"])


class TemplateDerivationTests(unittest.TestCase):
    """One example construct is friendlier input than two pasted strings."""

    def setUp(self) -> None:
        self.template = UPSTREAM + INSERT + DOWNSTREAM

    def test_flanks_are_recovered_from_an_example_construct(self) -> None:
        derived, matched = flanks.from_template(self.template, [INSERT], anchor_length=20)
        self.assertEqual(matched, INSERT)
        self.assertEqual(derived.upstream, UPSTREAM)
        self.assertEqual(derived.downstream, DOWNSTREAM)

    def test_the_right_insert_is_picked_out_of_many(self) -> None:
        others = ["TTTTTTTTTTTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGGGGGGGGGGG"]
        derived, matched = flanks.from_template(self.template, [*others, INSERT])
        self.assertEqual(matched, INSERT)
        self.assertEqual(derived.upstream, UPSTREAM)

    def test_a_reverse_complement_template_works(self) -> None:
        from nanopore3.sequence import reverse_complement

        derived, _ = flanks.from_template(reverse_complement(self.template), [INSERT])
        self.assertEqual(derived.upstream, UPSTREAM)

    def test_an_insert_with_a_few_errors_still_anchors_the_flanks(self) -> None:
        mutated = INSERT[:20] + "T" + INSERT[21:]
        derived, _ = flanks.from_template(self.template, [mutated])
        self.assertEqual(derived.upstream, UPSTREAM)

    def test_an_unrelated_template_is_refused_rather_than_guessed(self) -> None:
        with self.assertRaises(flanks.FlankError) as raised:
            flanks.from_template("ACGT" * 60, [INSERT])
        self.assertIn("identity", str(raised.exception))

    def test_a_template_that_starts_at_the_insert_is_refused(self) -> None:
        with self.assertRaises(flanks.FlankError) as raised:
            flanks.from_template(INSERT + DOWNSTREAM, [INSERT])
        self.assertIn("constant region", str(raised.exception))


class ExtractionRoundTripTests(unittest.TestCase):
    """The built reference must equal what the pipeline extracts from a read."""

    def test_a_synthetic_read_extracts_to_exactly_the_built_reference(self) -> None:
        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        built = f.flank(INSERT)
        read = "ACGTACGTACGT" + UPSTREAM + INSERT + DOWNSTREAM + "TGCATGCATGCA"
        extraction = extract_insert(read, f.left_anchor, f.right_anchor, max_edits=3)
        self.assertEqual(extraction.status, "found")
        self.assertEqual(extraction.sequence, built)

    def test_it_still_round_trips_when_the_read_carries_errors(self) -> None:
        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        read = "ACGT" + UPSTREAM[:5] + UPSTREAM[6:] + INSERT + DOWNSTREAM + "TGCA"
        extraction = extract_insert(read, f.left_anchor, f.right_anchor, max_edits=3)
        self.assertEqual(extraction.status, "found")
        self.assertIn(INSERT, extraction.sequence)

    def test_a_reverse_complement_read_extracts_the_forward_reference(self) -> None:
        from nanopore3.sequence import reverse_complement

        f = flanks.from_sequences(UPSTREAM, DOWNSTREAM, anchor_length=20)
        read = reverse_complement("ACGT" + UPSTREAM + INSERT + DOWNSTREAM + "TGCA")
        extraction = extract_insert(read, f.left_anchor, f.right_anchor, max_edits=3)
        self.assertEqual(extraction.sequence, f.flank(INSERT))


if __name__ == "__main__":
    unittest.main()
