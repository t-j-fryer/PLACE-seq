"""Positional reference matching: naming both parents of a chimeric clone."""

from __future__ import annotations

import random
import unittest

from nanopore3.assignment import ReferenceIndex
from nanopore3.chimera import (
    ChimeraGroup,
    ReadSignature,
    collapse_signature,
    group_chimeras,
    positional_profile,
    signature_for,
    synthesise_reference,
)


def make_references(count: int, length: int, seed: int) -> dict[str, str]:
    rng = random.Random(seed)
    return {
        f"design_{i}": "".join(rng.choice("ACGT") for _ in range(length))
        for i in range(count)
    }


def noisy(sequence: str, rate: float, rng: random.Random) -> str:
    out = []
    for base in sequence:
        roll = rng.random()
        if roll < rate / 3:
            continue
        if roll < 2 * rate / 3:
            out.append(rng.choice("ACGT"))
            continue
        out.append(base)
    return "".join(out)


class CollapseTests(unittest.TestCase):
    def test_consecutive_repeats_become_one_segment(self) -> None:
        self.assertEqual(collapse_signature(["a", "a", "a"]), ("a",))

    def test_unmatched_windows_are_skipped_not_treated_as_a_segment(self) -> None:
        """A window with no hit is missing evidence, not a third parent."""

        self.assertEqual(collapse_signature(["a", None, "a"]), ("a",))
        self.assertEqual(collapse_signature(["a", None, "b"]), ("a", "b"))

    def test_a_return_to_the_first_parent_is_kept_as_its_own_segment(self) -> None:
        self.assertEqual(collapse_signature(["a", "b", "a"]), ("a", "b", "a"))


class ProfileTests(unittest.TestCase):
    def setUp(self) -> None:
        self.refs = make_references(12, 600, seed=4)
        self.index = ReferenceIndex(self.refs, k=15, max_kmer_owners=10)

    def test_a_clean_read_gives_one_segment(self) -> None:
        rng = random.Random(1)
        read = noisy(self.refs["design_3"], 0.06, rng)
        result = signature_for("r", read, self.index)
        self.assertEqual(result.signature, ("design_3",))
        self.assertFalse(result.is_chimeric)

    def test_a_chimera_names_both_parents_in_order(self) -> None:
        rng = random.Random(2)
        left, right = self.refs["design_1"], self.refs["design_8"]
        read = noisy(left[:300] + right[300:], 0.06, rng)
        result = signature_for("r", read, self.index)
        self.assertEqual(result.signature, ("design_1", "design_8"))
        self.assertTrue(result.is_chimeric)
        self.assertTrue(result.junction_windows)

    def test_an_uneven_chimera_is_still_detected(self) -> None:
        """The case whole-read assignment gets confidently wrong."""

        rng = random.Random(3)
        left, right = self.refs["design_2"], self.refs["design_9"]
        read = noisy(left[:520] + right[520:], 0.06, rng)
        self.assertEqual(signature_for("r", read, self.index).signature,
                         ("design_2", "design_9"))

    def test_window_and_step_must_be_positive(self) -> None:
        with self.assertRaises(ValueError):
            positional_profile("ACGT", self.index, window=0)


class GroupingTests(unittest.TestCase):
    def _sig(self, read_id, signature):
        return ReadSignature(read_id, tuple(signature), tuple(signature))

    def test_reads_sharing_a_signature_form_one_clone(self) -> None:
        sigs = [self._sig(f"r{i}", ("a", "b")) for i in range(8)]
        groups, outcome = group_chimeras(sigs, minimum_depth=6)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0].size, 8)
        self.assertEqual(outcome["chimeric_read"], 8)

    def test_a_signature_below_the_depth_floor_is_noise_not_a_clone(self) -> None:
        """One noisy window can invent a segment; a real clone recurs."""

        sigs = [self._sig(f"r{i}", ("a", "b")) for i in range(3)]
        groups, outcome = group_chimeras(sigs, minimum_depth=6)
        self.assertEqual(groups, ())
        self.assertEqual(outcome["below_depth"], 3)

    def test_uniform_reads_are_counted_but_never_grouped(self) -> None:
        sigs = [self._sig(f"r{i}", ("a",)) for i in range(10)]
        groups, outcome = group_chimeras(sigs, minimum_depth=6)
        self.assertEqual(groups, ())
        self.assertEqual(outcome["uniform"], 10)

    def test_different_signatures_stay_separate_clones(self) -> None:
        sigs = ([self._sig(f"x{i}", ("a", "b")) for i in range(7)]
                + [self._sig(f"y{i}", ("a", "c")) for i in range(7)])
        groups, _ = group_chimeras(sigs, minimum_depth=6)
        self.assertEqual({g.signature for g in groups}, {("a", "b"), ("a", "c")})


class ScaffoldTests(unittest.TestCase):
    def test_two_parents_splice_at_the_junction(self) -> None:
        refs = {"a": "A" * 400, "b": "C" * 400}
        scaffold = synthesise_reference(("a", "b"), refs, junction_window=4, step=45)
        self.assertTrue(scaffold.startswith("A"))
        self.assertTrue(scaffold.endswith("C"))
        self.assertEqual(len(scaffold), 400)

    def test_more_than_two_parents_declines_to_splice(self) -> None:
        """Three junctions are not located well enough to trust a scaffold."""

        refs = {"a": "A" * 100, "b": "C" * 100, "c": "G" * 100}
        self.assertEqual(synthesise_reference(("a", "b", "c"), refs, 2), "")

    def test_an_unknown_parent_declines_to_splice(self) -> None:
        self.assertEqual(synthesise_reference(("a", "zz"), {"a": "A" * 50}, 1), "")

    def test_a_missing_junction_declines_to_splice(self) -> None:
        self.assertEqual(synthesise_reference(("a", "b"), {"a": "A"*50, "b": "C"*50}, None), "")


if __name__ == "__main__":
    unittest.main()
