from __future__ import annotations

import csv
import gzip
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from nanopore3.config import BarcodeSettings, load_config
from nanopore3.pipeline import _trimmed_barcodes, run_pipeline
from nanopore3.provenance import (
    StageDirectory,
    StageValidationError,
    canonical_digest,
    validate_stage_directory,
)

ROOT = Path(__file__).resolve().parents[1]


class ProvenancePipelineTests(unittest.TestCase):
    def test_barcode_trimming_matches_legacy_two_ended_rule(self) -> None:
        settings = BarcodeSettings(sequences={"x": "AACCGGTT"}, trim_bases=2)
        self.assertEqual(_trimmed_barcodes(settings), {"x": "CCGG"})

    def test_stage_publication_resume_and_corruption_detection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fingerprint = canonical_digest({"stage": "test"})
            with StageDirectory(root, "test", fingerprint) as stage:
                stage.output_path("result.txt").write_text("safe\n", encoding="utf-8")
            final = root / "stages" / "test"
            validate_stage_directory(final, expected_fingerprint=fingerprint)
            with StageDirectory(root, "test", fingerprint, resume=True) as stage:
                self.assertTrue(stage.reused)
            (final / "result.txt").write_text("changed\n", encoding="utf-8")
            with self.assertRaises(StageValidationError):
                validate_stage_directory(final)

    def test_end_to_end_synthetic_run_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            config = load_config(ROOT / "configs" / "example.yaml")
            run = run_pipeline(config, output_root=Path(directory), run_id="golden")
            self.assertTrue((run / "stages" / "06_report" / "report.html").is_file())
            assignment_summary = json.loads(
                (run / "stages" / "03_assignment" / "summary.json").read_text()
            )
            self.assertEqual(
                assignment_summary,
                {"assigned_alias_set": 3, "assigned_unique": 3, "no_match": 1},
            )
            with gzip.open(
                run / "stages" / "04_consensus" / "consensus.csv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(
                {row["reference_ids"] for row in rows},
                {"ref_A|ref_A_alias", "ref_B"},
            )
            with gzip.open(
                run / "stages" / "05_qc" / "qc.csv.gz",
                "rt",
                encoding="utf-8",
                newline="",
            ) as handle:
                qc_reader = csv.DictReader(handle)
                self.assertIn("alignment_identity", qc_reader.fieldnames or [])
                self.assertTrue(any(row["alignment_identity"] for row in qc_reader))
            resumed = run_pipeline(
                config, output_root=Path(directory), run_id="golden", resume=True
            )
            self.assertEqual(resumed, run)

    def test_parallel_worker_count_does_not_change_scientific_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            threaded = load_config(ROOT / "configs" / "example.yaml")
            serial = replace(
                threaded,
                parallel=replace(threaded.parallel, jobs=1, backend="serial"),
            )
            serial_run = run_pipeline(serial, output_root=root, run_id="serial")
            threaded_run = run_pipeline(threaded, output_root=root, run_id="threaded")
            artifacts = (
                "02_demux/demux_calls.csv.gz",
                "02_demux/demuxed_reads.jsonl.gz",
                "03_assignment/assignment_calls.csv.gz",
                "03_assignment/consensus_eligible.jsonl.gz",
                "04_consensus/consensus.csv.gz",
                "04_consensus/contributors.csv.gz",
                "04_consensus/consensus.fasta",
                "05_qc/qc.csv.gz",
            )
            for relative in artifacts:
                left = serial_run / "stages" / relative
                right = threaded_run / "stages" / relative
                if relative.endswith(".gz"):
                    with (
        gzip.open(left, "rb") as left_handle,
        gzip.open(right, "rb") as right_handle,
    ):
                        self.assertEqual(left_handle.read(), right_handle.read(), relative)
                else:
                    self.assertEqual(left.read_bytes(), right.read_bytes(), relative)


if __name__ == "__main__":
    unittest.main()
