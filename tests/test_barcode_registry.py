from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from nanopore3.barcodes import BarcodeRegistryError, read_barcode_panel
from nanopore3.config import ConfigError, load_config


ROOT = Path(__file__).resolve().parents[1]


class BarcodeRegistryTests(unittest.TestCase):
    def test_canonical_panels_have_expected_members_and_preserve_legacy_b1(self) -> None:
        reverse = read_barcode_panel(
            ROOT / "configs" / "barcodes" / "reverse_primer_families.csv",
            "rp_amplicon",
        )
        wells = read_barcode_panel(
            ROOT / "configs" / "barcodes" / "well_forward_primers_a1_h12.csv",
            "well_forward_full",
        )
        self.assertEqual(tuple(reverse.sequences), tuple(f"RP{i:02d}" for i in range(1, 13)))
        self.assertEqual(len(wells.sequences), 96)
        self.assertEqual(len(wells.sequences["B1"]), 59)
        self.assertEqual(len(wells.sequences["A1"]), 60)

    def test_config_resolves_registry_relative_to_yaml_and_records_digest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "panel.csv").write_text(
                "family_id,barcode_id,sequence\nplate,P1,AACCGGTT\n",
                encoding="utf-8",
            )
            (root / "run.yaml").write_text(
                "schema_version: 1\nrun_name: registry\noutput_root: runs\n"
                "inputs: [reads.fastq]\nreferences:\n  fasta: refs.fasta\n"
                "barcodes:\n  plate:\n    registry_csv: panel.csv\n"
                "    family_id: plate\n",
                encoding="utf-8",
            )
            config = load_config(root / "run.yaml")
            self.assertEqual(config.plate_barcodes.sequences, {"P1": "AACCGGTT"})
            self.assertEqual(len(config.plate_barcodes.registry_sha256 or ""), 64)
            self.assertEqual(config.plate_barcodes.registry_csv, (root / "panel.csv").resolve())

    def test_missing_family_and_mixed_inline_registry_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            panel = root / "panel.csv"
            panel.write_text(
                "family_id,barcode_id,sequence\nplate,P1,AACCGGTT\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(BarcodeRegistryError, "was not found"):
                read_barcode_panel(panel, "missing")
            config_path = root / "run.yaml"
            config_path.write_text(
                "schema_version: 1\nrun_name: bad\noutput_root: runs\n"
                "inputs: [reads.fastq]\nreferences:\n  fasta: refs.fasta\n"
                "barcodes:\n  plate:\n    sequences: {P1: AACCGGTT}\n"
                "    registry_csv: panel.csv\n    family_id: plate\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ConfigError, "either sequences or registry"):
                load_config(config_path)


if __name__ == "__main__":
    unittest.main()
