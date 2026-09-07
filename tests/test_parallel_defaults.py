"""Maximum-worker defaults and the single-CPU process fallback."""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from nanopore3.config import ParallelSettings, _parse_parallel, load_config
from nanopore3.pipeline import resolved_jobs
from nanopore3.runtime import plan_resources


class ParallelDefaultsTests(unittest.TestCase):
    def test_omitted_parallel_settings_use_all_cpus_and_processes(self):
        for settings in (ParallelSettings(), *map(_parse_parallel, (None, {}, {"jobs": 0}))):
            self.assertEqual(settings.backend, "process")
            self.assertEqual(settings.jobs, 0)
            self.assertEqual(settings.threads_per_job, 1)

    def test_shipped_example_resolves_to_all_detected_cpus(self):
        source = Path(__file__).resolve().parents[1] / "configs/example.yaml"
        config = load_config(source)
        with (patch("nanopore3.runtime.effective_cpu_count", return_value=16),
              patch("nanopore3.pipeline.effective_cpu_count", return_value=16)):
            self.assertEqual(config.parallel.backend, "process")
            plan = plan_resources(resolved_jobs(config), config.parallel.threads_per_job)
            self.assertEqual(plan.jobs, 16)

    def test_one_cpu_process_fallback_matches_serial_in_fresh_interpreters(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run([
                sys.executable, "-m", "nanopore3", "init", str(root / "example"),
            ], check=True, capture_output=True)
            source = root / "example/configs/example.yaml"
            body = yaml.safe_load(source.read_text(encoding="utf-8"))
            outputs = []
            for backend, jobs in (("serial", 1), ("process", 1), ("process", 0), ("process", 8)):
                name = f"{backend}-{jobs}"
                body["parallel"] = {"backend": backend, "jobs": jobs, "chunk_reads": 1}
                config_path = source.with_name(f"{name}.yaml")
                config_path.write_text(yaml.safe_dump(body), encoding="utf-8")
                # Patch CPU discovery inside a NEW interpreter: a preceding
                # initialized worker must not hide missing worker-local state.
                result = subprocess.run([
                    sys.executable, "-c",
                    "import os,sys; os.cpu_count=lambda:1; "
                    "from nanopore3.cli import main; raise SystemExit(main(sys.argv[1:]))",
                    "run", "--config", str(config_path), "--output", str(root / "runs"),
                    "--run-id", name, "--quiet",
                ], capture_output=True, text=True, timeout=60,
                    env={**os.environ, "MPLCONFIGDIR": str(root / "matplotlib")})
                self.assertEqual(result.returncode, 0, result.stderr)
                run = root / "runs" / name
                artifacts = []
                for relative in (
                    "02_demux/demux_calls.csv.gz", "03_assignment/assignment_calls.csv.gz",
                    "04_consensus/consensus.fasta", "05_qc/qc.csv.gz",
                ):
                    path = run / "stages" / relative
                    opener = gzip.open if path.suffix == ".gz" else open
                    with opener(path, "rb") as handle:
                        artifacts.append(handle.read())
                outputs.append(artifacts)
            self.assertTrue(all(artifacts == outputs[0] for artifacts in outputs))


if __name__ == "__main__":
    unittest.main()
