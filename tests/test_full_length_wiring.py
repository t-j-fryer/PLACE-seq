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


@unittest.skipUnless(FULL_LENGTH_CONFIG.is_file(), "shipped config not present")
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

    @unittest.skipUnless(INSERT_CONFIG.is_file(), "shipped config not present")
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
