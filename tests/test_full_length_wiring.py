"""Config-to-pipeline wiring for full-length references."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from nanopore3 import pipeline
from nanopore3.assignment import extract_insert
from nanopore3.config import ConfigError, ReferenceSettings, load_config
from nanopore3.flanks import from_sequences
from nanopore3.references import read_fasta, read_reference_libraries

UPSTREAM = "CATAATCCGCACGCATCTGGTCGATCCCGCGAAATTAATACGACTCACTATAGGGAGAGCAGCTATGCAGCTT"
DOWNSTREAM = (
    "AGTGGATCCGAAACACCAGGCACGTCGGAAAGCGCTACCCCTGAATCAGTGTTTACCTTGGAGGACTTTGTG"
    "CATCACCACCACCATCACTAATGACTCGAGTCTGGTAAAGAAACCGCTGCTGCGAAATTTGGATTGGCGAATGGGACGC"
)
INSERTS = {
    "d1": "ATGAAACGCGTGGAGGAAAAGGTCAAAGAGATTTTCGAAAAG",
    "d2": "ATGAAACGCGTGGAGGAAAAGGTCAAAGAGATTTTCGAATTT",
}


def write_fasta(directory: Path, name: str, records: dict[str, str]) -> Path:
    path = directory / name
    path.write_text(
        "".join(f">{key}\n{value}\n" for key, value in records.items()), encoding="ascii"
    )
    return path


class ConfigValidationTests(unittest.TestCase):
    def base(self, **kwargs) -> ReferenceSettings:
        return ReferenceSettings(fasta=(Path("x.fasta"),), **kwargs)

    def test_a_library_without_flanks_is_not_full_length(self) -> None:
        self.assertFalse(self.base().full_length)

    def test_explicit_flanks_make_a_library_full_length(self) -> None:
        settings = self.base(flanks_upstream=UPSTREAM, flanks_downstream=DOWNSTREAM)
        self.assertTrue(settings.full_length)

    def test_one_flank_alone_is_refused(self) -> None:
        with self.assertRaises(ConfigError):
            self.base(flanks_upstream=UPSTREAM)

    def test_flanks_and_a_template_together_are_refused(self) -> None:
        """Two definitions of the same constant regions can disagree."""

        with self.assertRaises(ConfigError):
            self.base(
                flanks_upstream=UPSTREAM,
                flanks_downstream=DOWNSTREAM,
                flanks_template=Path("t.fasta"),
            )

    def test_a_short_anchor_is_refused(self) -> None:
        with self.assertRaises(ConfigError):
            self.base(
                flanks_upstream=UPSTREAM, flanks_downstream=DOWNSTREAM, flanks_anchor_length=4
            )


class TransformTests(unittest.TestCase):
    """References are flanked as they are read, before digesting."""

    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = write_fasta(Path(self.directory.name), "inserts.fasta", INSERTS)
        self.flanks = from_sequences(UPSTREAM, DOWNSTREAM)

    def test_loaded_sequences_are_full_length(self) -> None:
        bundle = read_fasta(
            self.path, transform=lambda _id, sequence: self.flanks.flank(sequence)
        )
        for record in bundle.records:
            self.assertEqual(record.sequence, self.flanks.flank(INSERTS[record.id]))

    def test_the_digest_describes_what_is_searched_not_the_file(self) -> None:
        plain = read_fasta(self.path)
        flanked = read_fasta(
            self.path, transform=lambda _id, sequence: self.flanks.flank(sequence)
        )
        self.assertNotEqual(plain.digest, flanked.digest)

    def test_flanking_cannot_merge_two_distinct_inserts(self) -> None:
        """Alias groups must mean the same thing after flanking."""

        flanked = read_fasta(
            self.path, transform=lambda _id, sequence: self.flanks.flank(sequence)
        )
        self.assertEqual(flanked.alias_groups, ())

    def test_identical_inserts_stay_aliases_after_flanking(self) -> None:
        path = write_fasta(
            Path(self.directory.name), "dupes.fasta", {"a": INSERTS["d1"], "b": INSERTS["d1"]}
        )
        flanked = read_fasta(path, transform=lambda _id, s: self.flanks.flank(s))
        self.assertEqual(flanked.alias_groups, (("a", "b"),))

    def test_transforms_are_applied_per_library(self) -> None:
        collection = read_reference_libraries(
            {"flanked": self.path, "plain": self.path},
            transforms={"flanked": lambda _id, s: self.flanks.flank(s)},
        )
        flanked = collection.get("flanked").records[0]
        plain = collection.get("plain").records[0]
        self.assertGreater(len(flanked.sequence), len(plain.sequence))


class RegionSpanTests(unittest.TestCase):
    def test_spans_tile_the_reference_without_gaps_or_overlap(self) -> None:
        flanks = from_sequences(UPSTREAM, DOWNSTREAM)
        reference = flanks.flank(INSERTS["d1"])
        covered = sorted(pipeline.qc_regions(flanks, len(reference)).values())
        self.assertEqual(covered[0][0], 0)
        self.assertEqual(covered[-1][1], len(reference))
        for (_, end), (start, _) in zip(covered, covered[1:], strict=False):
            self.assertEqual(end, start)

    def test_the_insert_span_selects_the_insert(self) -> None:
        flanks = from_sequences(UPSTREAM, DOWNSTREAM)
        reference = flanks.flank(INSERTS["d2"])
        start, end = pipeline.qc_regions(flanks, len(reference))["insert"]
        self.assertEqual(reference[start:end], INSERTS["d2"])


FULL_LENGTH_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "runs"
    / "260608_full_length.yaml"
)
INSERT_CONFIG = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "runs"
    / "260608_sumo_lab_bs_fs_ph.yaml"
)


def shipped_data_available(config_path: Path) -> bool:
    """Whether a shipped config *and the references it points at* are present.

    The configs are tracked; the reference FASTAs they name are not - they live
    under ``runs/``, which is ignored, because they are experiment data rather
    than source.  Checking only the config passed locally and failed in CI, where
    the config exists and the FASTAs do not.
    """

    if not config_path.is_file():
        return False
    try:
        config = load_config(config_path)
    except Exception:  # an unloadable config is reported by ConfigValidationTests
        return False
    return all(
        path.is_file()
        for settings in config.reference_sets.values()
        for path in settings.fasta
    )


FULL_LENGTH_DATA = shipped_data_available(FULL_LENGTH_CONFIG)
INSERT_DATA = shipped_data_available(INSERT_CONFIG)


@unittest.skipUnless(FULL_LENGTH_DATA, "shipped references not present")
class ShippedConfigTests(unittest.TestCase):
    """The wiring is exercised through the configs the repository ships."""

    def test_flanks_resolve_for_every_declared_library(self) -> None:
        config = load_config(FULL_LENGTH_CONFIG)
        flanks = pipeline.resolve_flanks(config)
        self.assertEqual(sorted(flanks), ["lab", "sumo_ab"])
        for flank in flanks.values():
            self.assertTrue(flank.upstream.endswith("ATGCAGCTT"))
            self.assertGreater(flank.constant_bases, 800)

    def test_the_derived_anchors_become_the_region_motifs(self) -> None:
        config = load_config(FULL_LENGTH_CONFIG)
        flanks = pipeline.resolve_flanks(config)
        applied = pipeline.apply_flanks(config, flanks)
        for library_id, flank in flanks.items():
            settings = applied.reference_sets[library_id]
            self.assertEqual(settings.forward_motif, flank.left_anchor)
            self.assertEqual(settings.reverse_motif, flank.right_anchor)

    def test_references_load_at_full_length_through_the_pipeline_path(self) -> None:
        config = load_config(FULL_LENGTH_CONFIG)
        flanks = pipeline.resolve_flanks(config)
        collection = read_reference_libraries(
            {k: s.fasta for k, s in config.reference_sets.items()},
            transforms=pipeline.flank_transforms(flanks),
        )
        bundle = collection.get("sumo_ab")
        self.assertEqual(len(bundle.records), 684)
        shortest = min(len(record.sequence) for record in bundle.records)
        self.assertGreater(shortest, 900)

    def test_the_orf_spine_is_a_whole_number_of_codons(self) -> None:
        """A start-to-stop distance off a codon boundary would fail every clone."""

        config = load_config(FULL_LENGTH_CONFIG)
        spine = len(config.qc.upstream_constant) + len(config.qc.downstream_constant)
        self.assertEqual(spine % 3, 0)

    @unittest.skipUnless(INSERT_DATA, "shipped references not present")
    def test_an_insert_mode_config_resolves_no_flanks_and_is_unchanged(self) -> None:
        config = load_config(INSERT_CONFIG)
        self.assertEqual(pipeline.resolve_flanks(config), {})
        self.assertIs(pipeline.apply_flanks(config, {}), config)


if __name__ == "__main__":
    unittest.main()


class InsertViewTests(unittest.TestCase):
    """Chimera detection needs the discriminative part, not the whole amplicon."""

    def setUp(self) -> None:
        self.flanks = from_sequences(UPSTREAM, DOWNSTREAM)
        self.full = {
            "lib": {name: self.flanks.flank(seq) for name, seq in INSERTS.items()}
        }

    def test_the_insert_is_sliced_back_out_of_a_full_length_reference(self) -> None:
        references, _motifs = pipeline.insert_view(self.full, {"lib": self.flanks})
        self.assertEqual(references["lib"], INSERTS)

    def test_the_motifs_bound_the_insert_not_the_amplicon(self) -> None:
        _references, motifs = pipeline.insert_view(self.full, {"lib": self.flanks})
        left, right = motifs["lib"]
        self.assertTrue(UPSTREAM.endswith(left))
        self.assertTrue(DOWNSTREAM.startswith(right))
        self.assertNotEqual(left, self.flanks.left_anchor)

    def test_extracting_with_those_motifs_returns_the_insert(self) -> None:
        _references, motifs = pipeline.insert_view(self.full, {"lib": self.flanks})
        left, right = motifs["lib"]
        read = "ACGT" + UPSTREAM + INSERTS["d1"] + DOWNSTREAM + "TGCA"
        extraction = extract_insert(read, left, right, max_edits=3)
        self.assertEqual(extraction.sequence, INSERTS["d1"])

    def test_a_library_without_flanks_passes_through_untouched(self) -> None:
        references, motifs = pipeline.insert_view({"plain": INSERTS}, {})
        self.assertEqual(references["plain"], INSERTS)
        self.assertEqual(motifs, {})


class InsertGateTests(unittest.TestCase):
    """Acceptance floors belong where the references differ: the insert."""

    def setUp(self) -> None:
        self.flanks = from_sequences(UPSTREAM, DOWNSTREAM)
        inner_left, inner_right = self.flanks.insert_anchors()
        self.gate = pipeline.InsertGate(
            left_anchor=self.flanks.left_anchor,
            right_anchor=self.flanks.right_anchor,
            inner_left_anchor=inner_left,
            inner_right_anchor=inner_right,
            head=len(self.flanks.inner_upstream),
            tail=len(self.flanks.inner_downstream),
            motif_max_edits=3,
            minimum_identity=0.80,
            minimum_query_coverage=0.70,
            inserts=dict(INSERTS),
        )

    def read_of(self, insert: str) -> str:
        return "ACGT" + UPSTREAM + insert + DOWNSTREAM + "TGCA"

    def test_a_matching_clone_passes(self) -> None:
        identity, coverage, passed = self.gate.verdict(self.read_of(INSERTS["d1"]), "d1")
        self.assertTrue(passed)
        self.assertEqual(identity, 1.0)
        self.assertEqual(coverage, 1.0)

    def test_a_read_carrying_a_foreign_insert_segment_is_held_back(self) -> None:
        """The case that produced a spurious 100%-identity designed clone.

        The read is the assigned design preceded by a similar length of another
        design, so the amplicon still looks fine but the insert is half foreign.
        """

        foreign = "TTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGG"
        identity, coverage, passed = self.gate.verdict(
            self.read_of(foreign + INSERTS["d1"]), "d1"
        )
        self.assertFalse(passed)
        self.assertLess(coverage, 0.70)

    def test_the_insert_is_found_by_offset_when_its_boundary_is_damaged(self) -> None:
        """Assignable but not profilable reads must still be gated."""

        damaged = UPSTREAM[:-6] + "AAAAAA"  # inner anchor destroyed
        read = "ACGT" + damaged + INSERTS["d1"] + DOWNSTREAM + "TGCA"
        insert = self.gate.insert_of(read)
        self.assertIsNotNone(insert)
        self.assertIn(INSERTS["d1"][10:30], insert)

    def test_an_unusable_read_is_not_judged(self) -> None:
        self.assertIsNone(self.gate.verdict("ACGT" * 20, "d1"))

    def test_an_unknown_reference_is_not_judged(self) -> None:
        self.assertIsNone(self.gate.verdict(self.read_of(INSERTS["d1"]), "nope"))


@unittest.skipUnless(FULL_LENGTH_DATA, "shipped references not present")
class ShippedGateTests(unittest.TestCase):
    def test_gates_are_built_for_both_full_length_libraries(self) -> None:
        config = load_config(FULL_LENGTH_CONFIG)
        flanks = pipeline.resolve_flanks(config)
        config = pipeline.apply_flanks(config, flanks)
        collection = read_reference_libraries(
            {k: s.fasta for k, s in config.reference_sets.items()},
            transforms=pipeline.flank_transforms(flanks),
        )
        references = {
            library_id: {r.id: r.sequence for r in bundle.records}
            for library_id, bundle in collection.libraries
        }
        gates = pipeline.build_insert_gates(config, references)
        self.assertEqual(sorted(gates), ["lab", "sumo_ab"])
        gate = gates["sumo_ab"]
        self.assertEqual(len(gate.inserts), 684)
        # the floors are the library's own, applied where they discriminate
        self.assertEqual(gate.minimum_identity, config.reference_sets["sumo_ab"].minimum_identity)

    def test_a_library_can_opt_out(self) -> None:
        from dataclasses import replace as dc_replace

        config = load_config(FULL_LENGTH_CONFIG)
        flanks = pipeline.resolve_flanks(config)
        libraries = dict(config.reference_libraries)
        libraries["lab"] = dc_replace(libraries["lab"], insert_thresholds=False)
        config = dc_replace(config, reference_libraries=libraries)
        config = pipeline.apply_flanks(config, flanks)
        collection = read_reference_libraries(
            {k: s.fasta for k, s in config.reference_sets.items()},
            transforms=pipeline.flank_transforms(flanks),
        )
        references = {
            library_id: {r.id: r.sequence for r in bundle.records}
            for library_id, bundle in collection.libraries
        }
        self.assertEqual(sorted(pipeline.build_insert_gates(config, references)), ["sumo_ab"])
