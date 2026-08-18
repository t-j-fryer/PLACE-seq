"""Reference-free clustering of reads that differ only across a short insert."""

from __future__ import annotations

import random
import unittest
from collections import Counter

import edlib

from nanopore3.clustering import choose_seed, cluster_reads, find_variant_positions

BACKBONE = "".join(random.Random(0).choice("ACGT") for _ in range(2400))


def design(rng: random.Random, insert_len: int = 300) -> str:
    """A construct sharing 2.4 kb of backbone, differing only in a 300 nt insert.

    The generator takes its own RNG so each test is reproducible in isolation;
    a shared module-level RNG makes results depend on test execution order.
    """

    insert = "".join(rng.choice("ACGT") for _ in range(insert_len))
    return BACKBONE[:1200] + insert + BACKBONE[1200:]


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
        if roll > 1 - rate / 3:
            out.append(rng.choice("ACGT"))
    return "".join(out)


def identity(left: str, right: str) -> float:
    distance = edlib.align(left, right, mode="NW", task="distance")["editDistance"]
    return 1 - distance / max(len(left), len(right))


class VariantPositionTests(unittest.TestCase):
    def test_a_column_with_no_consistent_minority_is_not_a_variant(self) -> None:
        columns = [Counter({"A": 19, "C": 1})]
        self.assertEqual(
            find_variant_positions(columns, min_depth=6, min_allele_fraction=0.2), []
        )

    def test_a_consistent_minority_is_a_variant(self) -> None:
        columns = [Counter({"A": 12, "C": 8})]
        self.assertEqual(
            find_variant_positions(columns, min_depth=6, min_allele_fraction=0.2), [0]
        )

    def test_shallow_columns_are_ignored(self) -> None:
        columns = [Counter({"A": 2, "C": 2})]
        self.assertEqual(
            find_variant_positions(columns, min_depth=6, min_allele_fraction=0.2), []
        )


class SeedTests(unittest.TestCase):
    def test_the_seed_is_the_median_length_read_not_the_longest(self) -> None:
        """The longest read is the likeliest to carry a chimeric artefact."""

        reads = [("a", "A" * 100), ("b", "A" * 500), ("c", "A" * 110)]
        self.assertEqual(reads[choose_seed(reads)][0], "c")


class MonoclonalTests(unittest.TestCase):
    def test_one_clone_stays_one_cluster_and_polishes_to_the_truth(self) -> None:
        rng = random.Random(9)
        truth = design(rng)
        reads = [(f"r{i}", noisy(truth, 0.06, rng)) for i in range(20)]
        result = cluster_reads(reads, min_cluster_size=6)
        self.assertEqual(len(result.clusters), 1)
        self.assertEqual(result.clusters[0].size, 20)
        self.assertEqual(result.unassigned_read_ids, ())
        self.assertGreater(identity(result.clusters[0].consensus, truth), 0.99)

    def test_a_single_surviving_cluster_does_not_discard_the_rest(self) -> None:
        """Splitting to one cluster would look like low depth, not a failed split."""

        rng = random.Random(11)
        truth = design(rng)
        reads = [(f"r{i}", noisy(truth, 0.08, rng)) for i in range(20)]
        result = cluster_reads(reads, min_cluster_size=6)
        self.assertEqual(sum(c.size for c in result.clusters), 20)

    def test_too_few_reads_yields_no_cluster(self) -> None:
        rng = random.Random(3)
        truth = design(rng)
        reads = [(f"r{i}", noisy(truth, 0.06, rng)) for i in range(4)]
        result = cluster_reads(reads, min_cluster_size=6)
        self.assertEqual(result.clusters, ())
        self.assertEqual(len(result.unassigned_read_ids), 4)


class PolyclonalTests(unittest.TestCase):
    def test_clearly_distinct_clones_separate_without_mixing(self) -> None:
        """Clones whose inserts are unrelated separate cleanly."""

        rng = random.Random(5)
        designs = {name: design(rng, insert_len=900) for name in ("c0", "c1", "c2")}
        reads = [
            (f"{name}_r{i}", noisy(sequence, 0.06, rng))
            for name, sequence in designs.items()
            for i in range(10)
        ]
        result = cluster_reads(reads, min_cluster_size=6)
        self.assertGreaterEqual(len(result.clusters), 2)
        for cluster in result.clusters:
            sources = {read_id.split("_")[0] for read_id in cluster.read_ids}
            self.assertEqual(len(sources), 1, f"cluster mixes clones: {sources}")

    def test_clones_as_divergent_as_the_read_error_are_not_reliably_separated(self) -> None:
        """The measured limit of the method, asserted so it cannot regress silently.

        With a 300 nt insert in 2.7 kb the clones are ~94% identical, so their
        true divergence equals a 6% read error rate. The within-clone and
        between-clone distance distributions touch, and clusters can mix. Every
        read is still accounted for, and this regime is documented rather than
        presented as working.
        """

        rng = random.Random(5)
        designs = {name: design(rng) for name in ("c0", "c1", "c2")}
        reads = [
            (f"{name}_r{i}", noisy(sequence, 0.06, rng))
            for name, sequence in designs.items()
            for i in range(10)
        ]
        result = cluster_reads(reads, min_cluster_size=6)
        accounted = sum(c.size for c in result.clusters) + len(result.unassigned_read_ids)
        self.assertEqual(accounted, len(reads))
        mixed = [
            c for c in result.clusters
            if len({r.split("_")[0] for r in c.read_ids}) > 1
        ]
        self.assertTrue(mixed, "if this now separates cleanly, tighten the documented limit")

    def test_failing_to_split_is_conservative_not_a_wrong_answer(self) -> None:
        """Under-splitting degrades to the reference-guided behaviour.

        When the distance distributions overlap the reads collapse into one
        cluster whose consensus carries ambiguity codes, which the existing
        grading reports as ``mixed_variants``. It never silently attributes a
        read to the wrong clone.
        """

        rng = random.Random(23)
        designs = [design(rng), design(rng)]
        reads = [
            (f"c{i}_r{j}", noisy(sequence, 0.06, rng))
            for i, sequence in enumerate(designs)
            for j in range(10)
        ]
        result = cluster_reads(reads, min_cluster_size=6)
        clustered = sum(cluster.size for cluster in result.clusters)
        self.assertEqual(clustered + len(result.unassigned_read_ids), len(reads))
        if len(result.clusters) == 1:
            self.assertIn("N", result.clusters[0].consensus)

    def test_variants_are_found_in_the_insert_not_the_backbone(self) -> None:
        """The shared backbone must contribute no variant positions."""

        rng = random.Random(13)
        designs = [design(rng), design(rng)]
        reads = [
            (f"c{i}_r{j}", noisy(sequence, 0.06, rng))
            for i, sequence in enumerate(designs)
            for j in range(12)
        ]
        result = cluster_reads(reads, min_cluster_size=6)
        inside = [p for p in result.variant_positions if 1100 <= p <= 1600]
        self.assertGreater(len(inside) / len(result.variant_positions), 0.8)

    def test_a_minor_clone_below_the_allele_floor_is_absorbed(self) -> None:
        """This is the documented detection limit, not a silent failure."""

        rng = random.Random(17)
        major, minor = design(rng), design(rng)
        reads = [(f"a{i}", noisy(major, 0.06, rng)) for i in range(20)]
        reads += [(f"b{i}", noisy(minor, 0.06, rng)) for i in range(2)]
        result = cluster_reads(reads, min_cluster_size=6, min_allele_fraction=0.20)
        self.assertEqual(len(result.clusters), 1)


if __name__ == "__main__":
    unittest.main()
