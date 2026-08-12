"""Assignment must not depend on hash seeds, worker counts, or k-mer cascades."""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
import textwrap
import unittest

from nanopore3.assignment import ReferenceIndex, assign_read


def _synthetic_references(count: int, length: int, seed: int) -> dict[str, str]:
    rng = random.Random(seed)
    # A shared prefix and suffix give every reference k-mers that several
    # references own, which is what makes the specificity weighting meaningful.
    shared_head = "".join(rng.choice("ACGT") for _ in range(40))
    shared_tail = "".join(rng.choice("ACGT") for _ in range(40))
    return {
        f"REF{index:03d}": shared_head
        + "".join(rng.choice("ACGT") for _ in range(length))
        + shared_tail
        for index in range(count)
    }


def _mutate(sequence: str, rate: float, seed: int) -> str:
    rng = random.Random(seed)
    bases = []
    for base in sequence:
        roll = rng.random()
        if roll < rate / 3:
            continue
        if roll < 2 * rate / 3:
            bases.append(rng.choice("ACGT"))
            continue
        bases.append(base)
        if roll > 1 - rate / 3:
            bases.append(rng.choice("ACGT"))
    return "".join(bases)


_HASH_SEED_PROBE = textwrap.dedent(
    """
    import json, sys
    sys.path.insert(0, {src!r})
    from tests.test_assignment_determinism import _mutate, _synthetic_references
    from nanopore3.assignment import ReferenceIndex, assign_read

    references = _synthetic_references(60, 300, seed=11)
    indexes = tuple(
        ReferenceIndex(references, k=k, max_kmer_owners=10) for k in (15, 11, 9)
    )
    results = []
    for read_index, reference_id in enumerate(sorted(references)):
        query = _mutate(references[reference_id], 0.06, seed=read_index)
        # Scores are compared as exact hexadecimal floats.  Rounding them would
        # hide the last-place differences that reorder candidates and move them
        # across the minimum-k-mer-score cutoff.
        shortlists = [
            [
                [candidate.aliases[0], candidate.kmer_score.hex(), candidate.matching_kmers]
                for candidate in index.shortlist(query, top_n=10, min_kmer_score=0.0)
            ]
            for index in indexes
        ]
        call, k = assign_read(
            query, indexes, top_n=5, min_kmer_score=2.0,
            min_identity=0.80, min_query_coverage=0.70,
            min_reference_coverage=0.70, min_identity_margin=0.02,
        )
        results.append([
            call.status,
            list(call.reference_ids),
            k,
            None if call.best is None else call.best.kmer_score.hex(),
            None if call.identity_margin is None else call.identity_margin.hex(),
            shortlists,
        ])
    print(json.dumps(results))
    """
)


class AssignmentDeterminismTests(unittest.TestCase):
    def setUp(self) -> None:
        self.references = _synthetic_references(60, 300, seed=11)
        self.indexes = tuple(
            ReferenceIndex(self.references, k=k, max_kmer_owners=10) for k in (15, 11, 9)
        )

    def test_shortlist_scores_do_not_depend_on_kmer_iteration_order(self) -> None:
        """Scores are exact sums over specificity classes, not float accumulations."""

        index = self.indexes[0]
        query = _mutate(self.references["REF007"], 0.06, seed=3)
        expected = index.shortlist(query, top_n=8, min_kmer_score=0.0)
        # Re-instantiating with the references in a different insertion order
        # changes internal dict and set layout without changing the biology.
        shuffled_ids = sorted(self.references, reverse=True)
        shuffled = ReferenceIndex(
            {reference_id: self.references[reference_id] for reference_id in shuffled_ids},
            k=index.k,
            max_kmer_owners=10,
        )
        self.assertEqual(shuffled.shortlist(query, top_n=8, min_kmer_score=0.0), expected)

    def test_assignment_is_stable_across_hash_seeds(self) -> None:
        """Process workers inherit different PYTHONHASHSEED values than the parent."""

        source_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        script = _HASH_SEED_PROBE.format(src=source_root)
        outputs = []
        for hash_seed in ("0", "1", "12345"):
            environment = dict(os.environ, PYTHONHASHSEED=hash_seed)
            completed = subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True,
                text=True,
                check=True,
                env=environment,
                cwd=source_root,
            )
            outputs.append(json.loads(completed.stdout))
        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(outputs[0], outputs[2])
        self.assertTrue(any(row[0] == "assigned_unique" for row in outputs[0]))

    def test_kmer_rescue_matches_exhaustive_alignment(self) -> None:
        """The bounded k-mer rescue must not lose calls the full sweep would make."""

        thresholds = dict(
            top_n=5,
            min_kmer_score=2.0,
            min_identity=0.80,
            min_query_coverage=0.70,
            min_reference_coverage=0.70,
            min_identity_margin=0.02,
        )
        assigned = 0
        for read_index, reference_id in enumerate(sorted(self.references)):
            query = _mutate(self.references[reference_id], 0.06, seed=read_index)
            kmer_call, _ = assign_read(query, self.indexes, rescue="kmer", **thresholds)
            all_call, _ = assign_read(query, self.indexes, rescue="all", **thresholds)
            self.assertEqual(kmer_call.status, all_call.status)
            self.assertEqual(kmer_call.reference_ids, all_call.reference_ids)
            if kmer_call.status == "assigned_unique":
                assigned += 1
                self.assertEqual(kmer_call.reference_ids, (reference_id,))
        self.assertGreater(assigned, len(self.references) // 2)

    def test_cascade_reports_the_k_that_produced_the_call(self) -> None:
        query = _mutate(self.references["REF042"], 0.05, seed=99)
        call, k = assign_read(
            query,
            self.indexes,
            top_n=5,
            min_kmer_score=2.0,
            min_identity=0.80,
            min_query_coverage=0.70,
            min_reference_coverage=0.70,
        )
        self.assertIn(k, {index.k for index in self.indexes})
        self.assertEqual(call.status, "assigned_unique")
        self.assertEqual(call.reference_ids, ("REF042",))

    def test_descending_kmer_sizes_are_required_by_the_cascade(self) -> None:
        with self.assertRaises(ValueError):
            assign_read("ACGT" * 30, ())


if __name__ == "__main__":
    unittest.main()
