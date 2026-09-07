"""Recomputing the analysis stages from a finished run's intermediates.

The property that matters is not that reuse is fast but that it is refused when it
would be wrong: a stage is inherited only if the current configuration produces the
fingerprint that stage recorded.  These tests pin the refusals, not just the
happy path.
"""

from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from nanopore3.provenance import analysis_implementation
from nanopore3.rerun import (
    INPUT_CONSUMING_STAGES,
    STAGE_ORDER,
    RerunError,
    completed_stages,
    prepare_rerun,
    stages_before,
)


def make_run(root: Path, stages: tuple[str, ...], *, run_id: str = "source") -> Path:
    run = root / run_id
    (run / "stages").mkdir(parents=True)
    (run / "run.json").write_text(
        json.dumps(
            {
                "run_id": run_id,
                "analysis_implementation": analysis_implementation(),
                "preflight": {
                    "inputs": [{"sample_id": "s", "sha256": "abc", "path": "/gone/reads.fastq"}]
                },
            }
        ),
        encoding="utf-8",
    )
    for stage in stages:
        directory = run / "stages" / stage
        directory.mkdir()
        payload = stage.encode("ascii")
        (directory / "data.txt").write_bytes(payload)
        (directory / "manifest.json").write_text(
            json.dumps(
                {
                    "schema_version": 1,
                    "stage": stage,
                    "fingerprint": f"fp-{stage}",
                    "pipeline_version": "0.3.0",
                    "created_utc": "2026-08-21T00:00:00Z",
                    "completed_utc": "2026-08-21T00:00:01Z",
                    "parameters": {"schema": 1},
                    "input_digests": {"s": "abc"},
                    "backend_versions": {},
                    "runtime": {"analysis_implementation": analysis_implementation()},
                    "artifacts": [
                        {
                            "path": "data.txt",
                            "sha256": hashlib.sha256(payload).hexdigest(),
                            "size_bytes": len(payload),
                            "media_type": "text/plain",
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        (directory / "_SUCCESS").write_text(f"fp-{stage}\n", encoding="utf-8")
    return run


class StageOrderTests(unittest.TestCase):
    def test_stages_before_is_everything_earlier(self) -> None:
        self.assertEqual(
            stages_before("04_consensus"),
            ("01_ingest", "02_demux", "03_assignment", "03b_chimera"),
        )

    def test_the_first_stage_has_nothing_before_it(self) -> None:
        self.assertEqual(stages_before(STAGE_ORDER[0]), ())

    def test_an_unknown_stage_lists_the_real_ones(self) -> None:
        with self.assertRaises(RerunError) as caught:
            stages_before("04_consensuss")
        self.assertIn("04_consensus", str(caught.exception))

    def test_the_input_consuming_stages_come_first(self) -> None:
        # What makes a detached rerun possible: inherit both and the FASTQ is
        # never opened.
        self.assertEqual(STAGE_ORDER[: len(INPUT_CONSUMING_STAGES)], INPUT_CONSUMING_STAGES)


class PrepareRerunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def test_completed_stages_are_reported_in_pipeline_order(self) -> None:
        run = make_run(self.root, ("02_demux", "01_ingest"))
        self.assertEqual(completed_stages(run), ("01_ingest", "02_demux"))

    def test_earlier_stages_are_carried_into_a_new_directory(self) -> None:
        source = make_run(self.root, ("01_ingest", "02_demux", "03_assignment"))
        destination = self.root / "new"
        record = prepare_rerun(source, destination, "04_consensus")
        self.assertEqual(
            record["inherited_stages"], ["01_ingest", "02_demux", "03_assignment"]
        )
        self.assertEqual(record["recomputed_from"], "04_consensus")
        self.assertTrue((destination / "stages" / "02_demux" / "data.txt").is_file())
        # The source is never modified: a rerun writes a new run.
        self.assertTrue((source / "stages" / "02_demux" / "data.txt").is_file())

    def test_artifacts_are_hard_linked_rather_than_duplicated(self) -> None:
        source = make_run(self.root, ("01_ingest", "02_demux"))
        destination = self.root / "new"
        prepare_rerun(source, destination, "03_assignment")
        original = (source / "stages" / "02_demux" / "data.txt").stat()
        carried = (destination / "stages" / "02_demux" / "data.txt").stat()
        self.assertEqual(original.st_ino, carried.st_ino)

    def test_the_input_digest_is_carried_forward_without_the_file(self) -> None:
        # The point of the feature: provenance still names the FASTQ and its
        # checksum even though nothing opened it.
        source = make_run(self.root, ("01_ingest", "02_demux"))
        record = prepare_rerun(source, self.root / "new", "03_assignment")
        inputs = record["source_preflight"]["inputs"]
        self.assertEqual(inputs[0]["sha256"], "abc")

    def test_an_optional_stage_that_never_ran_is_not_an_error(self) -> None:
        # Chimera detection only runs when enabled, so its absence is normal.
        source = make_run(self.root, ("01_ingest", "02_demux", "03_assignment"))
        record = prepare_rerun(source, self.root / "new", "04_consensus")
        self.assertEqual(record["not_inherited"], ["03b_chimera"])
        self.assertNotIn("03b_chimera", record["inherited_stages"])

    def test_a_fingerprint_mismatch_is_refused_by_name(self) -> None:
        source = make_run(self.root, ("01_ingest", "02_demux"))
        with self.assertRaises(RerunError) as caught:
            prepare_rerun(
                source,
                self.root / "new",
                "03_assignment",
                expected_fingerprints={"02_demux": "fp-different"},
            )
        self.assertIn("02_demux", str(caught.exception))

    def test_a_source_with_nothing_usable_is_refused(self) -> None:
        source = make_run(self.root, ("06_report",))
        with self.assertRaises(RerunError):
            prepare_rerun(source, self.root / "new", "04_consensus")

    def test_a_directory_that_is_not_a_run_is_refused(self) -> None:
        (self.root / "notarun").mkdir()
        with self.assertRaises(RerunError) as caught:
            prepare_rerun(self.root / "notarun", self.root / "new", "04_consensus")
        self.assertIn("run.json", str(caught.exception))

    def test_rerunning_from_the_first_stage_is_a_normal_run(self) -> None:
        source = make_run(self.root, ("01_ingest",))
        with self.assertRaises(RerunError) as caught:
            prepare_rerun(source, self.root / "new", "01_ingest")
        self.assertIn("nothing to inherit", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
