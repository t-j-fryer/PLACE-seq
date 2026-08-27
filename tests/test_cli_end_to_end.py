"""Every command line surface, exercised the way a user reaches it.

The unit tests in this suite cover functions. Three defects in a row lived in the
paths *between* functions - a second call site resolving a value differently, a
measurement that was absent rather than zero, an exception class the entry point did
not catch - and every one of them had a passing unit test. So each command is run
here as a subprocess, and the two things asserted of every failure are the two that
kept being wrong: a non-zero exit, and a sentence rather than a traceback.
"""

from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
_STATE: dict[str, Path] = {}


def run_cli(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", "nanopore3", *args],
        capture_output=True,
        text=True,
        cwd=str(cwd or REPO),
    )


def setUpModule() -> None:
    tmp = tempfile.TemporaryDirectory()
    _STATE["tmp"] = tmp  # type: ignore[assignment]
    root = Path(tmp.name)
    _STATE["root"] = root
    assert run_cli("init", str(root / "ex")).returncode == 0
    _STATE["config"] = root / "ex" / "configs" / "example.yaml"
    assert run_cli(
        "run", "--config", str(_STATE["config"]), "--output", str(root / "runs"),
        "--run-id", "base", "--quiet",
    ).returncode == 0
    _STATE["base"] = root / "runs" / "base"


def tearDownModule() -> None:
    _STATE["tmp"].cleanup()  # type: ignore[attr-defined]


class CleanFailureTests(unittest.TestCase):
    """No user mistake may produce a traceback, whatever the command."""

    def variant(self, name: str, mutate) -> Path:
        body = yaml.safe_load(_STATE["config"].read_text(encoding="utf-8"))
        mutate(body)
        path = _STATE["config"].parent / f"{name}.yaml"
        path.write_text(yaml.safe_dump(body), encoding="utf-8")
        return path

    def assert_clean_failure(self, result, hint: str = "") -> None:
        combined = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, f"expected failure: {combined[:200]}")
        self.assertNotIn("Traceback", combined, combined[-400:])
        self.assertIn("error:", combined, combined[-400:])
        if hint:
            self.assertIn(hint, combined)

    def test_an_unknown_configuration_key(self) -> None:
        path = self.variant("unknown_key", lambda b: b.__setitem__("bogus", 1))
        self.assert_clean_failure(run_cli("validate", "--config", str(path)), "bogus")

    def test_an_out_of_range_value(self) -> None:
        path = self.variant(
            "bad_identity", lambda b: b["references"].__setitem__("minimum_identity", 1.8)
        )
        self.assert_clean_failure(run_cli("validate", "--config", str(path)), "between")

    def test_a_missing_input_file(self) -> None:
        path = self.variant(
            "missing_input",
            lambda b: b["inputs"][0].__setitem__("path", "../fixtures/absent.fastq"),
        )
        self.assert_clean_failure(run_cli("validate", "--config", str(path)), "missing")

    def test_an_unset_path_variable(self) -> None:
        path = self.variant(
            "unset_var",
            lambda b: b["inputs"][0].__setitem__("path", "${NOT_SET_ANYWHERE}/r.fastq"),
        )
        self.assert_clean_failure(
            run_cli("validate", "--config", str(path)), "NOT_SET_ANYWHERE"
        )

    def test_a_consensus_backend_whose_tools_are_absent(self) -> None:
        # Reached the user as a traceback until Nanopore3Error became the base of
        # every user-facing error. A very likely first-run mistake.
        path = self.variant(
            "no_backend", lambda b: b["consensus"].__setitem__("backend", "mafft_spoa")
        )
        result = run_cli("validate", "--config", str(path))
        if result.returncode == 0:
            self.skipTest("MAFFT and SPOA are installed on this machine")
        self.assert_clean_failure(result, "mafft_spoa")

    def test_a_missing_configuration_file(self) -> None:
        self.assert_clean_failure(run_cli("validate", "--config", "/nonexistent.yaml"))

    def test_an_unknown_rerun_stage(self) -> None:
        # Also reached the user as a traceback.
        self.assert_clean_failure(
            run_cli(
                "rerun", "--config", str(_STATE["config"]),
                "--from-run", str(_STATE["base"]), "--from", "nonsense",
                "--run-id", "never",
            ),
            "04_consensus",
        )

    def test_a_run_directory_that_already_exists(self) -> None:
        self.assert_clean_failure(
            run_cli(
                "run", "--config", str(_STATE["config"]),
                "--output", str(_STATE["root"] / "runs"), "--run-id", "base",
            ),
            "already exists",
        )

    def test_subsample_on_something_that_is_not_a_fastq(self) -> None:
        notes = _STATE["root"] / "notes.txt"
        notes.write_text("not a fastq\n", encoding="ascii")
        self.assert_clean_failure(
            run_cli(
                "subsample", "--input", str(notes),
                "--output", str(_STATE["root"] / "out.fastq"), "--reads", "1",
            ),
            "not a FASTQ",
        )

    def test_a_layout_with_no_reference_map(self) -> None:
        def mutate(body):
            body["reference_libraries"] = {
                "a": {"fasta": "../fixtures/references.fasta"},
                "b": {"fasta": "../fixtures/references.fasta"},
            }
            body.pop("references", None)

        path = self.variant("two_libraries", mutate)
        self.assert_clean_failure(
            run_cli("validate", "--config", str(path)), "plate_reference_map"
        )


class SucceedingCommandTests(unittest.TestCase):
    def test_doctor_reports_and_emits_json(self) -> None:
        self.assertEqual(run_cli("doctor").returncode, 0)
        result = run_cli("doctor", "--json")
        self.assertEqual(result.returncode, 0)
        report = json.loads(result.stdout)
        self.assertIn("optional_binaries", report)

    def test_validate_reports_progress_and_can_be_silenced(self) -> None:
        loud = run_cli("validate", "--config", str(_STATE["config"]))
        self.assertEqual(loud.returncode, 0)
        self.assertIn("preflight", loud.stderr)
        quiet = run_cli("validate", "--config", str(_STATE["config"]), "--quiet")
        self.assertEqual(quiet.returncode, 0)
        self.assertNotIn("preflight", quiet.stderr)

    def test_quick_validation_skips_the_record_count(self) -> None:
        result = run_cli("validate", "--config", str(_STATE["config"]), "--quick")
        self.assertEqual(result.returncode, 0)
        self.assertNotIn("counting records", result.stderr)

    def test_a_run_names_its_three_outputs_and_records_provenance(self) -> None:
        base = _STATE["base"]
        for relative in (
            "run.json",
            "consensus_by_plate/index.csv",
            "stages/06_report/report.html",
        ):
            self.assertTrue((base / relative).is_file(), relative)
        metadata = json.loads((base / "run.json").read_text(encoding="utf-8"))
        for key in ("config", "config_digest", "preflight", "runtime"):
            self.assertIn(key, metadata)
        self.assertTrue(metadata["preflight"]["inputs"][0]["sha256"])

    def test_a_run_reports_each_stage(self) -> None:
        result = run_cli(
            "run", "--config", str(_STATE["config"]),
            "--output", str(_STATE["root"] / "runs"), "--run-id", "loud",
        )
        self.assertEqual(result.returncode, 0)
        for stage in ("01_ingest", "02_demux", "05_qc", "06_report"):
            self.assertIn(f"stage {stage}", result.stderr)

    def test_layout_reports_when_pooling_is_not_configured(self) -> None:
        result = run_cli("layout", "--config", str(_STATE["config"]))
        self.assertEqual(result.returncode, 0)
        self.assertIn("not enabled", result.stdout)

    def test_subsample_takes_the_head_of_a_file(self) -> None:
        out = _STATE["root"] / "sub.fastq.gz"
        result = run_cli(
            "subsample", "--input", str(_STATE["root"] / "ex/fixtures/reads.fastq"),
            "--output", str(out), "--reads", "2",
        )
        self.assertEqual(result.returncode, 0)
        with gzip.open(out, "rt", encoding="ascii") as handle:
            self.assertEqual(handle.read().count("@"), 2)

    def test_init_refuses_to_overwrite(self) -> None:
        result = run_cli("init", str(_STATE["root"] / "ex"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("error:", result.stdout + result.stderr)


class RerunEndToEndTests(unittest.TestCase):
    def test_it_works_from_every_stage_that_has_a_predecessor(self) -> None:
        for index, stage in enumerate(
            ("02_demux", "03_assignment", "04_consensus", "05_qc", "06_report")
        ):
            result = run_cli(
                "rerun", "--config", str(_STATE["config"]),
                "--from-run", str(_STATE["base"]),
                "--from", stage, "--output", str(_STATE["root"] / "runs"),
                "--run-id", f"rr{index}", "--quiet",
            )
            self.assertEqual(result.returncode, 0, f"{stage}: {result.stderr[-200:]}")
            self.assertTrue(
                (_STATE["root"] / "runs" / f"rr{index}" / "run.json").is_file()
            )

    def test_the_first_stage_has_nothing_to_inherit(self) -> None:
        result = run_cli(
            "rerun", "--config", str(_STATE["config"]),
            "--from-run", str(_STATE["base"]), "--from", "01_ingest",
            "--run-id", "never2",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("nothing to inherit", result.stdout + result.stderr)

    def test_a_change_reaching_an_inherited_stage_is_refused_and_cleans_up(self) -> None:
        body = yaml.safe_load(_STATE["config"].read_text(encoding="utf-8"))
        body["library"]["minimum_read_length"] = int(
            body["library"]["minimum_read_length"]
        ) + 5
        changed = _STATE["config"].parent / "demux_changed.yaml"
        changed.write_text(yaml.safe_dump(body), encoding="utf-8")

        destination = _STATE["root"] / "runs" / "refused"
        result = run_cli(
            "rerun", "--config", str(changed), "--from-run", str(_STATE["base"]),
            "--from", "04_consensus", "--output", str(_STATE["root"] / "runs"),
            "--run-id", "refused", "--quiet",
        )
        combined = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, combined[:200])
        self.assertNotIn("Traceback", combined)
        self.assertIn("earliest stage your change affects", combined)
        # A refusal must not leave a half-built run behind.
        self.assertFalse(destination.exists())

    def test_it_runs_with_the_input_file_detached(self) -> None:
        # The reason the command exists: the FASTQ is not needed once the stages
        # that read it are inherited.
        fastq = _STATE["root"] / "ex/fixtures/reads.fastq"
        hidden = fastq.with_suffix(".hidden")
        fastq.rename(hidden)
        try:
            result = run_cli(
                "rerun", "--config", str(_STATE["config"]),
                "--from-run", str(_STATE["base"]), "--from", "04_consensus",
                "--output", str(_STATE["root"] / "runs"), "--run-id", "detached",
                "--quiet",
            )
            self.assertEqual(result.returncode, 0, result.stderr[-300:])
            # And a plain run of the same configuration must still fail.
            plain = run_cli(
                "run", "--config", str(_STATE["config"]),
                "--output", str(_STATE["root"] / "runs"), "--run-id", "plain_detached",
            )
            self.assertNotEqual(plain.returncode, 0)
        finally:
            hidden.rename(fastq)


class DeterminismTests(unittest.TestCase):
    def test_worker_count_and_backend_do_not_change_the_science(self) -> None:
        import hashlib

        digests = set()
        for index, (jobs, backend) in enumerate(
            ((1, "serial"), (2, "process"), (0, "process"), (2, "thread"))
        ):
            body = yaml.safe_load(_STATE["config"].read_text(encoding="utf-8"))
            body["parallel"] = {"jobs": jobs, "backend": backend, "chunk_reads": 1}
            path = _STATE["config"].parent / f"par{index}.yaml"
            path.write_text(yaml.safe_dump(body), encoding="utf-8")
            result = run_cli(
                "run", "--config", str(path), "--output", str(_STATE["root"] / "runs"),
                "--run-id", f"par{index}", "--quiet",
            )
            self.assertEqual(result.returncode, 0, f"{jobs}/{backend}: {result.stderr[-200:]}")
            run = _STATE["root"] / "runs" / f"par{index}"
            with gzip.open(run / "stages/05_qc/qc.csv.gz", "rb") as handle:
                qc = hashlib.sha256(handle.read()).hexdigest()
            fasta = hashlib.sha256(
                (run / "stages/04_consensus/consensus.fasta").read_bytes()
            ).hexdigest()
            digests.add((fasta, qc))
        self.assertEqual(len(digests), 1, "outputs differ by worker count or backend")


if __name__ == "__main__":
    unittest.main()
