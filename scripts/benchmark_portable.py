#!/usr/bin/env python3
"""Reproducible synthetic whole-pipeline timing, memory and output equivalence.

Install psutil separately (or use the MCP extra). No private inputs are needed.
This measures a small reference panel; it does not predict production throughput.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import random
import statistics
import subprocess
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import psutil
import yaml


def measure(command: list[str], log: Path) -> dict:
    started = time.perf_counter()
    peak = 0
    with log.open("wb") as handle:
        process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT)
        parent = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                memory = sum(
                    p.memory_info().rss for p in [parent, *parent.children(recursive=True)]
                )
                peak = max(peak, memory)
            except psutil.NoSuchProcess:
                pass
            time.sleep(0.02)
        code = process.wait()
    return {
        "seconds": time.perf_counter() - started,
        "peak_tree_rss_mib": peak / 2**20, "returncode": code,
    }


def digest(path: Path) -> str:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new benchmark directory")
    parser.add_argument("--reads", type=int, default=12_000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--process-jobs", type=int, nargs="+", default=[2, 4, 8, 0],
        help="process worker counts to sweep; 0 requests all detected CPUs",
    )
    args = parser.parse_args()
    if args.reads < 1 or args.repeats < 1:
        parser.error("reads and repeats must be positive")
    if any(jobs < 0 for jobs in args.process_jobs):
        parser.error("process-jobs must be nonnegative")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    subprocess.run([sys.executable, "-m", "nanopore3", "init", str(root / "example")], check=True)
    config_path = root / "example/configs/example.yaml"
    original = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    reads = root / "example/fixtures/reads.fastq"
    lines = reads.read_text(encoding="ascii").splitlines(keepends=True)
    records = [lines[i:i + 4] for i in range(0, len(lines), 4)]
    with reads.open("w", encoding="ascii", newline="\n") as handle:
        for index in range(args.reads):
            record = records[index % len(records)]
            handle.write(f"@synthetic-{index}\n" + "".join(record[1:]))
    results = []
    modes = [("serial", 1), ("thread", 4)] + [
        ("process", jobs) for jobs in dict.fromkeys(args.process_jobs)
    ]
    paths = {}
    for backend, jobs in modes:
        config = {**original, "parallel": {
            "backend": backend, "jobs": jobs, "threads_per_job": 1, "chunk_reads": 1000,
        }}
        path = config_path.with_name(f"{backend}-{jobs}.yaml")
        path.write_text(yaml.safe_dump(config), encoding="utf-8")
        paths[(backend, jobs)] = path
    order = random.Random(0)
    for repeat in range(args.repeats):
        sweep = modes.copy()
        order.shuffle(sweep)
        for backend, jobs in sweep:
            path = paths[(backend, jobs)]
            run_id = f"{backend}-{jobs}-{repeat}"
            result = measure([
                sys.executable, "-m", "nanopore3", "run", "--config", str(path),
                "--output", str(root / "runs"), "--run-id", run_id, "--quiet",
            ], root / f"{run_id}.log")
            run = root / "runs" / run_id
            result.update(backend=backend, requested_jobs=jobs, repeat=repeat,
                          reads_per_second=(args.reads / result["seconds"]
                                            if result["returncode"] == 0 else None))
            if (run / "run.json").is_file():
                metadata = json.loads((run / "run.json").read_text(encoding="utf-8"))
                result["resolved_jobs"] = metadata["preflight"]["resources"]["jobs"]
            if result["returncode"] == 0:
                result["digests"] = {relative: digest(run / relative) for relative in (
                    "stages/02_demux/demux_calls.csv.gz",
                    "stages/03_assignment/assignment_calls.csv.gz",
                    "stages/04_consensus/consensus.fasta",
                    "stages/05_qc/qc.csv.gz",
                )}
            else:
                result["error_tail"] = (root / f"{run_id}.log").read_text(
                    encoding="utf-8", errors="replace",
                )[-1500:]
            results.append(result)
            (root / "trials.json").write_text(json.dumps(results, indent=2) + "\n")
            print(json.dumps({k: v for k, v in result.items() if k != "digests"}), flush=True)
    succeeded = [row for row in results if row["returncode"] == 0]
    equivalent = bool(succeeded) and all(
        row["digests"] == succeeded[0]["digests"] for row in succeeded
    )
    all_succeeded = len(succeeded) == len(results)
    reporting = {}
    for name in ("matplotlib", "numpy", "pandas", "scipy"):
        try:
            reporting[name] = version(name)
        except PackageNotFoundError:
            reporting[name] = None
    report = {
        "schema_version": 2, "order_seed": 0,
        "platform": platform.platform(), "python": platform.python_version(),
        "cpu_count": psutil.cpu_count(), "reads": args.reads, "repeats": args.repeats,
        "packages": {name: version(name) for name in ("nanopore3", "edlib", "PyYAML", "psutil")},
        "reporting_packages": reporting,
        "memory_method": (
            "20 ms sampling; summed process-tree RSS, includes shared pages per process"
        ),
        "equivalent_successful_runs": equivalent, "all_runs_succeeded": all_succeeded,
        "results": results,
        "medians": {f"{backend}-{jobs}": {
            field: statistics.median(
                row[field] for row in succeeded
                if row["backend"] == backend and row["requested_jobs"] == jobs
            )
            for field in ("seconds", "reads_per_second", "peak_tree_rss_mib")
        } for backend, jobs in modes if any(
            row["backend"] == backend and row["requested_jobs"] == jobs for row in succeeded
        )},
    }
    (root / "results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Successful output equivalence: {equivalent}; all runs succeeded: {all_succeeded}")
    print(f"Results: {root / 'results.json'}")
    return 0 if equivalent and all_succeeded else 1


if __name__ == "__main__":
    raise SystemExit(main())
