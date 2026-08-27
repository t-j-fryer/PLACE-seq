from __future__ import annotations

import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from nanopore3.config import ConfigError, load_config
from nanopore3.pipeline import run_pipeline
from nanopore3.references import FastaFormatError, read_fasta, read_reference_libraries


class ReferenceLibraryConfigTests(unittest.TestCase):
    def _write(self, root: Path, relative: str, text: str) -> Path:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def _base_yaml(self, reference_section: str) -> str:
        return (
            "schema_version: 1\n"
            "run_name: reference-test\n"
            "output_root: results\n"
            "inputs:\n"
            "  - path: reads.fastq\n"
            f"{reference_section}"
        )

    def test_legacy_single_references_remains_available_as_default(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = self._write(
                root,
                "config.yaml",
                self._base_yaml("references:\n  fasta: refs.fasta\n"),
            )
            config = load_config(config_path)
            self.assertIsNotNone(config.references)
            self.assertEqual(tuple(config.reference_sets), ("default",))
            self.assertEqual(config.reference_library_id_for_plate("RP01"), "default")
            self.assertIs(
                config.reference_library_for_plate("RP07"), config.references
            )

    def test_named_libraries_and_plate_map_are_resolved_stably(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config_path = self._write(
                root,
                "config.yaml",
                self._base_yaml(
                    "reference_libraries:\n"
                    "  zeta:\n"
                    "    fasta: zeta.fasta\n"
                    "    kmer_sizes: [21, 15]\n"
                    "  alpha:\n"
                    "    fasta: [alpha-a.fasta, alpha-b.fasta]\n"
                    "plate_reference_map:\n"
                    "  RP07: zeta\n"
                    "  RP01: alpha\n"
                ),
            )
            config = load_config(config_path)
            self.assertIsNone(config.references)
            self.assertEqual(tuple(config.reference_sets), ("alpha", "zeta"))
            self.assertEqual(
                config.reference_library_for_plate("RP01").fasta,
                (
                    (root / "alpha-a.fasta").resolve(),
                    (root / "alpha-b.fasta").resolve(),
                ),
            )
            self.assertEqual(
                config.reference_library_for_plate("RP07").kmer_sizes, (21, 15)
            )
            with self.assertRaisesRegex(KeyError, "RP02"):
                config.reference_library_for_plate("RP02")
            serialized = config.as_dict()
            self.assertNotIn("source_path", serialized)
            self.assertEqual(
                tuple(serialized["reference_libraries"]), ("alpha", "zeta")
            )
            self.assertEqual(
                tuple(serialized["plate_reference_map"]), ("RP01", "RP07")
            )

    def test_map_target_and_duplicate_paths_are_validated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bad_target = self._write(
                root,
                "bad-target.yaml",
                self._base_yaml(
                    "reference_libraries:\n"
                    "  alpha:\n"
                    "    fasta: refs.fasta\n"
                    "plate_reference_map:\n"
                    "  RP01: missing\n"
                ),
            )
            with self.assertRaisesRegex(ConfigError, "unknown reference library"):
                load_config(bad_target)

            duplicate_path = self._write(
                root,
                "duplicate-path.yaml",
                self._base_yaml(
                    "reference_libraries:\n"
                    "  alpha:\n"
                    "    fasta: [refs.fasta, ./refs.fasta]\n"
                ),
            )
            with self.assertRaisesRegex(ConfigError, "duplicate paths"):
                load_config(duplicate_path)


class ReferenceLibraryReadTests(unittest.TestCase):
    def test_duplicate_ids_are_rejected_within_but_allowed_across_libraries(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.fasta"
            second = root / "second.fasta"
            first.write_text(">shared\nACGTACGT\n", encoding="utf-8")
            second.write_text(">shared\nTTTTCCCC\n", encoding="utf-8")

            with self.assertRaisesRegex(FastaFormatError, "Duplicate reference ID"):
                read_fasta((first, second))

            collection = read_reference_libraries(
                {"second": second, "first": first}
            )
            self.assertEqual(collection.ids, ("first", "second"))
            self.assertEqual(collection.get("first").ids, ("shared",))
            self.assertNotEqual(
                collection.get("first").digest,
                collection.get("second").digest,
            )

    def test_duplicate_paths_are_rejected_and_collection_digest_is_stable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.fasta"
            second = root / "second.fasta"
            first.write_text(">a\nACGTACGT\n", encoding="utf-8")
            second.write_text(">b\nTTTTCCCC\n", encoding="utf-8")
            with self.assertRaisesRegex(FastaFormatError, "must not occur more than once"):
                read_fasta((first, first))
            forward = read_reference_libraries({"a": first, "b": second})
            reverse = read_reference_libraries({"b": second, "a": first})
            self.assertEqual(forward.ids, reverse.ids)
            self.assertEqual(forward.digest, reverse.digest)


class ReferenceLibraryPipelineTests(unittest.TestCase):
    def test_assignment_is_routed_by_plate_and_unmapped_plate_is_explicit(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "alpha.fasta").write_text(
                ">shared\nATGAAACCCGGGTTTAAACCCGGGTTTTAA\n", encoding="utf-8"
            )
            (root / "beta.fasta").write_text(
                ">shared\nATGTTTGGGCCCAAATTTGGGCCCAAATAA\n", encoding="utf-8"
            )
            inserts = {
                "alpha": ("AAAACCCC", "ATGAAACCCGGGTTTAAACCCGGGTTTTAA"),
                "beta": ("AAGGCCAT", "ATGTTTGGGCCCAAATTTGGGCCCAAATAA"),
                "unmapped": ("TTGACCAA", "ATGAAACCCGGGTTTAAACCCGGGTTTTAA"),
            }
            fastq_parts: list[str] = []
            for read_id, (barcode, insert) in inserts.items():
                sequence = barcode + "GATTACA" + insert + "CCGGAAT"
                fastq_parts.append(
                    f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                )
            (root / "reads.fastq").write_text("".join(fastq_parts), encoding="ascii")
            (root / "run.yaml").write_text(
                "schema_version: 1\n"
                "run_name: routed\n"
                "output_root: runs\n"
                "inputs:\n  - path: reads.fastq\n    sample_id: routed\n"
                "reference_libraries:\n"
                "  alpha:\n    fasta: alpha.fasta\n    kmer_sizes: [5]\n"
                "  beta:\n    fasta: beta.fasta\n    kmer_sizes: [5]\n"
                "plate_reference_map:\n  P1: alpha\n  P2: beta\n"
                "library:\n  forward_motif: GATTACA\n  reverse_motif: CCGGAAT\n"
                "barcodes:\n  plate:\n    sequences:\n"
                "      P1: AAAACCCC\n      P2: AAGGCCAT\n      P3: TTGACCAA\n"
                "    max_edits: 0\n    search_window: 8\n"
                "    search_ends: [head]\n    minimum_margin: 1\n"
                "consensus:\n  minimum_depth: 1\n",
                encoding="utf-8",
            )
            config = load_config(root / "run.yaml")
            run = run_pipeline(config, output_root=root / "output", run_id="routed")
            with gzip.open(
                run / "stages" / "03_assignment" / "assignment_calls.csv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                rows = {row["original_read_id"]: row for row in csv.DictReader(handle)}
            self.assertEqual(rows["alpha"]["reference_library_id"], "alpha")
            self.assertEqual(rows["beta"]["reference_library_id"], "beta")
            self.assertEqual(rows["alpha"]["reference_ids"], "shared")
            self.assertEqual(rows["beta"]["reference_ids"], "shared")
            self.assertEqual(
                rows["unmapped"]["assignment_status"],
                "unmapped_reference_library",
            )


if __name__ == "__main__":
    unittest.main()
