"""Regression coverage for data loss, stale reuse, uncapped reads and resources."""

from __future__ import annotations

import csv
import gzip
import json
import os
import random
import sqlite3
import tracemalloc
import unittest
from collections import Counter
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from nanopore3.cli import _subsample, main
from nanopore3.config import ConfigError, _parse_parallel, load_config
from nanopore3.consensus import ConsensusRead, build_reference_consensus, select_reads
from nanopore3.grouping import DiskStore, GroupMemoryError, GroupStorageError
from nanopore3.pipeline import PipelineError, _detect_chimeras, run_pipeline
from nanopore3.provenance import (
    StageValidationError,
    analysis_implementation,
    compute_stage_fingerprint,
    validate_stage_directory,
)
from nanopore3.rerun import prepare_rerun
from nanopore3.runtime import effective_cpu_count, group_memory_limit, plan_resources

ROOT = Path(__file__).resolve().parents[1]
RECORD = b"@read description\r\nACGT\r\n+read description\r\nIIII\r\n"


class TemporaryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)


class SubsampleSafetyTests(TemporaryTest):
    def test_existing_output_and_all_input_aliases_are_preserved(self):
        source = self.root / "input.fastq"
        source.write_bytes(RECORD)
        outputs = [source, self.root / "existing.fastq", self.root / "hard.fastq"]
        outputs[1].write_bytes(b"existing data")
        os.link(source, outputs[2])
        for path in outputs:
            before = path.read_bytes()
            with self.assertRaises(PipelineError):
                _subsample(source, path, 1)
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(source.read_bytes(), RECORD)
        for name, target in (("symbolic", source), ("dangling", self.root / "missing")):
            alias = self.root / name
            try:
                alias.symlink_to(target)
            except OSError:
                continue  # Windows can disallow symlinks without developer mode.
            with self.assertRaises(PipelineError):
                _subsample(source, alias, 1)
            self.assertTrue(alias.is_symlink())

    def test_malformed_prefix_never_publishes_partial_output(self):
        source, dest = self.root / "in.fastq", self.root / "out.fastq.gz"
        for bad in (
            b"@bad\nAC\n+\nI\n",
            b"@bad\nAC\n-\nII\n",
            b"@bad\nAC\n+\n",
            b"@bad\nAC\n+\nI \n",
            b"@bad\nAZ\n+\nII\n",
        ):
            source.write_bytes(RECORD + bad)
            with self.assertRaises(PipelineError):
                _subsample(source, dest, 2)
            self.assertFalse(dest.exists())
            self.assertEqual(source.read_bytes(), RECORD + bad)
            self.assertFalse(list(self.root.glob(".subsample-*")))
        # Only the requested prefix is validated; raw bytes remain unchanged.
        self.assertEqual(_subsample(source, dest, 1), 1)
        self.assertEqual(gzip.decompress(dest.read_bytes()), RECORD)

    def test_racing_destination_and_unsupported_link_fail_without_clobber(self):
        source, dest = self.root / "in.fastq", self.root / "out.fastq"
        source.write_bytes(RECORD)
        original_link = os.link

        def race(temp, destination):
            destination.write_bytes(b"other writer")
            original_link(temp, destination)

        with patch("nanopore3.cli.os.link", side_effect=race):
            with self.assertRaises(FileExistsError):
                _subsample(source, dest, 1)
        self.assertEqual(dest.read_bytes(), b"other writer")
        dest.unlink()
        with patch("nanopore3.cli.os.link", side_effect=OSError("unsupported filesystem")):
            with self.assertRaises(OSError):
                _subsample(source, dest, 1)
        self.assertFalse(dest.exists())
        self.assertFalse(list(self.root.glob(".subsample-*")))

    def test_truncated_gzip_has_a_clean_cli_error_and_no_output(self):
        source, dest = self.root / "in.fastq.gz", self.root / "out.fastq"
        source.write_bytes(gzip.compress(RECORD)[:-8])
        self.assertEqual(
            main(["subsample", "--input", str(source), "--output", str(dest), "--reads", "2"]), 2
        )
        self.assertFalse(dest.exists())
        self.assertFalse(list(self.root.glob(".subsample-*")))


class ImplementationTests(TemporaryTest):
    def test_resident_interpreter_refuses_changed_source_before_preflight(self):
        config = load_config(ROOT / "configs/example.yaml")
        with (
            patch("nanopore3.provenance._IMPORTED_IMPLEMENTATION", {"sha256": "old"}),
            patch("nanopore3.pipeline.validate_inputs") as preflight,
        ):
            with self.assertRaisesRegex(StageValidationError, "restart Python/the MCP server"):
                run_pipeline(config, output_root=self.root, run_id="changed")
            preflight.assert_not_called()
        self.assertFalse((self.root / "changed").exists())

    def test_content_identity_is_portable_and_tracks_edits_without_git(self):
        a, b = self.root / "a", self.root / "b"
        for root, newline in ((a, b"\n"), (b, b"\r\n")):
            (root / "presets").mkdir(parents=True)
            (root / "core.py").write_bytes(b"value = 1" + newline)
            (root / "presets/default.yaml").write_bytes(b"jobs: 0" + newline)
        first = analysis_implementation(a)
        self.assertEqual(first, analysis_implementation(b))
        (a / "core.py").write_text("value = 2\n")
        changed = analysis_implementation(a)
        self.assertNotEqual(first, changed)
        with patch("nanopore3.provenance.analysis_implementation", return_value=first):
            fingerprint = compute_stage_fingerprint("stage", pipeline_version="same")
        with patch("nanopore3.provenance.analysis_implementation", return_value=changed):
            self.assertNotEqual(
                fingerprint, compute_stage_fingerprint("stage", pipeline_version="same")
            )
        (b / "README.md").write_text("Documentation only")
        self.assertEqual(first, analysis_implementation(b))

    def test_resume_and_rerun_refuse_legacy_or_changed_identity_before_work(self):
        config = load_config(ROOT / "configs/example.yaml")
        run = self.root / "source"
        run.mkdir()
        for identity in (None, {"scheme": "nanopore3-source-v1", "sha256": "old"}):
            metadata = {} if identity is None else {"analysis_implementation": identity}
            path = run / "run.json"
            path.write_text(json.dumps(metadata))
            before = path.read_bytes()
            with patch("nanopore3.pipeline.validate_inputs") as preflight:
                with self.assertRaisesRegex(StageValidationError, "resume/rerun refused"):
                    run_pipeline(config, output_root=self.root, run_id="source", resume=True)
                preflight.assert_not_called()
            with self.assertRaisesRegex(StageValidationError, "resume/rerun refused"):
                prepare_rerun(run, self.root / "new", "04_consensus")
            self.assertFalse((self.root / "new").exists())
            self.assertEqual(path.read_bytes(), before)


class ResourceTests(TemporaryTest):
    def cgroup(self, version=2):
        proc = self.root / "proc"
        (proc / "self").mkdir(parents=True)
        mount = self.root / "cgroup"
        child = mount / "job"
        child.mkdir(parents=True)
        if version == 2:
            membership = "0::/job\n"
            tail = "cgroup2 cgroup rw"
        else:
            membership = "3:cpu,cpuset,memory:/job\n"
            tail = "cgroup cgroup rw,cpu,cpuset,memory"
        (proc / "self/cgroup").write_text(membership)
        (proc / "self/mountinfo").write_text(f"10 9 0:1 / {mount} rw - {tail}\n")
        return proc, mount, child

    @patch("os.cpu_count", return_value=32)
    @patch("os.sched_getaffinity", return_value=set(range(12)), create=True)
    @patch("os.process_cpu_count", return_value=20, create=True)
    def test_v2_parent_quota_cpuset_and_nested_thread_budget(self, *_):
        proc, mount, child = self.cgroup()
        (child / "cpu.max").write_text("max 100000")
        (mount / "cpu.max").write_text("350000 100000")
        (child / "cpuset.cpus.effective").write_text("1-8,10-11")
        self.assertEqual(effective_cpu_count(proc=proc), 3)
        with patch("nanopore3.runtime.effective_cpu_count", return_value=3):
            plan = plan_resources(32, 64)
            self.assertEqual((plan.jobs, plan.threads_per_job, plan.available_cpus), (1, 3, 3))
            self.assertLessEqual(plan_resources(32, 2).jobs * 2, 3)
        (mount / "cpu.max").write_text("50000 100000")
        self.assertEqual(effective_cpu_count(proc=proc), 1)
        (mount / "memory.max").write_text(str(256 * 1024 * 1024))
        self.assertEqual(group_memory_limit(512, proc=proc), 64 * 1024 * 1024)

    @patch("os.cpu_count", return_value=32)
    @patch("os.sched_getaffinity", return_value=set(range(8)), create=True)
    @patch("os.process_cpu_count", return_value=32, create=True)
    def test_v1_quota_and_cpuset_and_missing_controls(self, *_):
        proc, mount, child = self.cgroup(version=1)
        (child / "cpu.cfs_quota_us").write_text("400000")
        (child / "cpu.cfs_period_us").write_text("100000")
        (mount / "cpuset.cpus").write_text("2,5-6")
        self.assertEqual(effective_cpu_count(proc=proc), 3)
        (mount / "memory.limit_in_bytes").write_text(str(128 * 1024 * 1024))
        self.assertEqual(group_memory_limit(512, proc=proc), 32 * 1024 * 1024)
        self.assertEqual(effective_cpu_count(proc=self.root / "missing"), 8)

    def test_memory_configuration_is_validated(self):
        self.assertEqual(_parse_parallel({}).group_memory_mb, 512)
        self.assertEqual(_parse_parallel({"group_memory_mb": 32}).group_memory_mb, 32)
        for bad in (0, -1, True, "lots"):
            with self.assertRaises(ConfigError):
                _parse_parallel({"group_memory_mb": bad})


class GroupingTests(TemporaryTest):
    def test_suspended_cursors_close_on_storage_failure(self):
        with self.assertRaisesRegex(GroupStorageError, "check free disk space"):
            with DiskStore(self.root) as store:
                for index in range(3):
                    store.append("reads", index, key=("group",))
                groups = store.groups("reads")
                _, rows = next(groups)
                self.assertEqual(next(rows), 0)
                raise sqlite3.OperationalError("database or disk is full")
        self.assertFalse(store._cursors)
        self.assertEqual(list(self.root.iterdir()), [])
        groups.close()
        rows.close()

    def test_disk_groups_keep_tuple_order_and_input_order(self):
        keys = [
            ("a", ("x",)),
            ("a", ("x", "y")),
            ('a"', ("x",)),
            ("a\\", ("x",)),
            ("é", ("x",)),
            ("a", ("x\x00",)),
        ]
        with DiskStore(self.root) as store:
            for i in range(3):
                for key in reversed(keys):
                    store.append("reads", i, key=key)
            actual = [(key, list(rows)) for key, rows in store.groups("reads")]
            self.assertEqual(actual, [(key, [0, 1, 2]) for key in sorted(keys)])
        self.assertEqual(list(self.root.iterdir()), [])

    def test_many_groups_stream_and_cleanup_on_failure(self):
        tracemalloc.start()
        try:
            with self.assertRaisesRegex(RuntimeError, "injected"):
                with DiskStore(self.root) as store:
                    for i in range(2000):
                        store.append("reads", "A" * 4096, key=(f"group-{i}",))
                    count = 0
                    for _key, rows in store.groups("reads"):
                        self.assertEqual(sum(len(row) for row in rows), 4096)
                        count += 1
                    self.assertEqual(count, 2000)
                    self.assertLess(tracemalloc.get_traced_memory()[1], 2 * 1024 * 1024)
                    raise RuntimeError("injected")
        finally:
            tracemalloc.stop()
        self.assertEqual(list(self.root.iterdir()), [])


class UnlimitedConsensusTests(TemporaryTest):
    def config(self):
        config = load_config(ROOT / "configs/example.yaml")
        return replace(
            config,
            parallel=replace(config.parallel, jobs=1, backend="serial"),
            consensus=replace(config.consensus, maximum_reads=0),
        )

    def test_public_zero_cap_means_all_and_remains_order_independent(self):
        reads = [ConsensusRead(f"r{i}", "ACGT") for i in range(5)]
        selected = select_reads(reads, max_reads=0, seed=1, group_id="g")
        self.assertEqual(set(r.read_uid for r in selected), {r.read_uid for r in reads})
        a = build_reference_consensus("ACGT", reads, group_id="g", max_reads=0, min_depth=1)
        b = build_reference_consensus(
            "ACGT", reads[::-1], group_id="g", max_reads=5, min_depth=1
        )
        self.assertEqual(a, b)
        self.assertEqual(a.n_reads_used, 5)

    def test_pipeline_uses_every_eligible_read_and_passes_bounded_threads(self):
        config = self.config()
        config = replace(config, parallel=replace(config.parallel, threads_per_job=999))
        with (
            patch(
                "nanopore3.pipeline.build_reference_consensus", wraps=build_reference_consensus
            ) as build,
            patch("nanopore3.runtime.effective_cpu_count", return_value=2),
        ):
            run = run_pipeline(config, output_root=self.root, run_id="all")
        self.assertTrue(build.call_count)
        self.assertTrue(all(call.kwargs["threads"] == 2 for call in build.call_args_list))
        with gzip.open(run / "stages/04_consensus/consensus.csv.gz", "rt") as handle:
            rows = list(csv.DictReader(handle))
        self.assertTrue(rows)
        self.assertEqual(
            json.loads((run / "stages/04_consensus/summary.json").read_text()),
            dict(Counter(row["status"] for row in rows)),
        )
        self.assertTrue(all(row["n_reads_used"] == row["n_reads_available"] for row in rows))
        self.assertEqual(sum(int(row["n_reads_used"]) for row in rows), 6)
        eligible = run / "stages/03_assignment/consensus_eligible.jsonl.gz"
        with gzip.open(eligible, "rt") as handle:
            ids = {json.loads(line)["read_uid"] for line in handle}
        with gzip.open(run / "stages/04_consensus/contributors.csv.gz", "rt") as handle:
            self.assertEqual({row["read_uid"] for row in csv.DictReader(handle)}, ids)
        manifest = validate_stage_directory(run / "stages/04_consensus")
        self.assertEqual(manifest.runtime["analysis_implementation"], analysis_implementation())
        self.assertFalse(list(run.rglob("*.sqlite")))
        with patch("nanopore3.pipeline.group_memory_limit", return_value=1):
            with self.assertRaisesRegex(GroupMemoryError, "No reads were silently dropped"):
                run_pipeline(config, output_root=self.root, run_id="limited")
        self.assertFalse((self.root / "limited/stages/04_consensus").exists())
        self.assertFalse(list((self.root / "limited").rglob("*.sqlite")))

    def test_chimera_zero_cap_writes_clone_and_all_claimed_reads(self):
        config = self.config()
        config = replace(
            config,
            chimera=replace(
                config.chimera,
                enabled=True,
                minimum_depth=2,
                write_pcr_origin=True,
                exclude_reads="all",
            ),
        )
        rng = random.Random(14)
        refs = {name: "".join(rng.choice("ACGT") for _ in range(600)) for name in ("a", "b")}
        sequence = refs["a"][:300] + refs["b"][300:]
        library = config.reference_library_id_for_plate("plate_1")
        demux = self.root / "reads.gz"
        with gzip.open(demux, "wt") as handle:
            for i in range(5):
                handle.write(
                    json.dumps(
                        {
                            "plate_id": "plate_1",
                            "well_id": "A1",
                            "read_uid": f"r{i}",
                            "sequence": sequence,
                        }
                    )
                    + "\n"
                )
        outputs = []
        for cap in (5, 0):
            output = self.root / f"cap-{cap}"
            with patch(
                "nanopore3.pipeline.build_reference_consensus", wraps=build_reference_consensus
            ) as build:
                summary = _detect_chimeras(
                    replace(config, consensus=replace(config.consensus, maximum_reads=cap)),
                    demux,
                    {library: refs},
                    output,
                    extractors={library: lambda s: s},
                )
            self.assertEqual(summary["written"], 1)
            self.assertEqual(summary["reads_claimed"], 5)
            self.assertEqual(len(build.call_args.args[1]), 5)
            outputs.append(
                {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}
            )
        self.assertEqual(outputs[0], outputs[1])
        with patch("nanopore3.pipeline.group_memory_limit", return_value=1):
            with self.assertRaises(GroupMemoryError):
                _detect_chimeras(
                    config,
                    demux,
                    {library: refs},
                    self.root / "limited",
                    extractors={library: lambda s: s},
                )
        self.assertFalse(list(self.root.rglob("*.sqlite")))
