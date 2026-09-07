# Repository assessment — 2026-09-07

**Subsequent change:** at the user's request, the package defaults, shipped examples
and run profiles now use `backend: process`, `jobs: 0`, `threads_per_job: 1`.
The one-resolved-worker assignment failure discovered below has been fixed and
covered by fresh-interpreter equivalence tests. The benchmark tables retain
their original pre-fix evidence; references below to eight-worker production
settings describe the configuration at the time of that assessment.

PLACE-seq (formerly Nanopore3) has a useful, well-tested portable core and unusually strong output
provenance for a pre-alpha research pipeline. It is suitable for controlled
analyses and iterative validation. It is not yet a resource-bounded workflow
engine or a fully locked, production-validated scientific application.

**Four-priority follow-up completed:** subsampling now protects existing files,
source identity controls reuse, a zero consensus cap includes all eligible reads,
and effective CPU allocation plus disk-backed grouping/estimated memory guards
replace the problematic resource behavior. See the updated findings and
[configuration guide](configuration.md). A broad rewrite or GPU port is not
supported by the current evidence.

**Interface follow-up:** CLI help/diagnostics, the local Jupyter/Colab notebook
and MCP now share the updated resource/reuse behavior. MCP has 13 tools including
safe subsampling. The interface suite passes 443 tests; installed-wheel stdio
verification passes. See the latest lab entry and [MCP guide](mcp.md).

**Hosted validation follow-up:** [GitHub Actions run 34153711406](https://github.com/t-j-fryer/PLACE-seq/actions/runs/34153711406)
passed Ubuntu/Python 3.10 and 3.13, macOS/Python 3.12, Windows/Python 3.12 and the
Linux installed-wheel/report/MCP workflow for implementation commit `827da1e`.
Earlier statements about local-only verification below describe their original
audit stage. Live Colab remains unverified. See the
[publication guide](updates/2026-09-07-portability-and-mcp.md) for the platform table.

## Evidence and scope

- Read the current package, configuration, CLI, stage/provenance implementation,
  tests, benchmark scripts, CI, notebook, architecture and recent lab entries.
- Baseline on macOS ARM64 / Python 3.12.10: **414 tests passed in 19.85 s**;
  `ruff check src scripts tests` passed before changes.
- Measured a fresh, reproducible synthetic workload, with three repeats each
  for serial, four threads and four processes. Preserved raw measurements and
  artifact hashes in [benchmark evidence](benchmarks/2026-09-07-synthetic.json).
- Reproduced the uncapped-consensus failure and CPU-budget violation described
  below. They are existing core defects, left open in this change.
- Added and exercised the MCP adapter over actual stdio and in-process HTTP
  protocol sessions, including complete synthetic execution and cancellation.
- Final local suite: **420 passed in 22.41 s**; lint, compilation and `pip check`
  passed. A real loopback HTTP server completed the SDK smoke workflow. A clean
  wheel environment outside the checkout ran the core example without MCP,
  then ran the stdio smoke workflow after installing the MCP extra.
- Linux, Windows, WSL and hosted Colab were **not executed in this session**.
  A configured CI matrix is coverage intent, not a successful test result.
  Cross-platform CI has been extended; its remote results remain to be observed.
- No new real sequencing run or scientific accuracy study was performed.
  Existing full-scale numbers in the lab notebook are historical evidence.

## Structure

| Area | Assessment |
|---|---|
| Packaging | Standard `src/` layout, wheel assets, CLI and Python API; just edlib/PyYAML in the core |
| Configuration | Strict unknown-key rejection, dataclasses, presets, named path variables, explicit plate routing |
| Analysis | Useful domain modules for I/O, matching, references, consensus, QC, flanks and pooled layouts |
| Orchestration | `pipeline.py` is 2,434 lines and mixes scheduling, I/O, analysis decisions and report assembly |
| Configuration implementation | `config.py` is 1,580 lines; parsing and semantic validation are tightly combined |
| Evidence | Checksummed atomic stages, per-read decisions, contributors and explicit QC uncertainty are strengths |
| Tests | Broad synthetic and regression coverage, subprocess CLI tests, spawn/backend determinism tests |
| Utilities | Useful diagnostics, but several import private pipeline functions and can drift from the real run path |
| Documentation | Worked examples and newest-first lab history are strong; some architecture promises describe a target, not current behavior |

The two largest modules account for about one third of the 12,073 lines in the
pre-MCP package. Extract stage execution functions behind explicit inputs/results
before adding more features to `run_pipeline`. Keep demultiplexing, assignment
and QC behavior under the existing golden tests during that extraction.

The MCP layer now lives separately in `mcp_server.py` (protocol and tools) and
`mcp_jobs.py` (workspace and subprocess lifecycle). It invokes the same CLI used
by humans and notebooks. Typed tools, resources and a workflow prompt make it
discoverable without putting protocol code in the scientific modules.

## Measured performance

Command from the repo, with the package and psutil installed:

```bash
python scripts/benchmark_portable.py \
  --output runs/portability-audit-20260907 --reads 12000 --repeats 3 --process-jobs 4
```

Use a new output directory for a repeat. The script creates synthetic inputs
from the packaged fixture, assigns unique read headers, and runs the entire CLI
including preflight, compressed outputs, QC and available report generation.
It compares decompressed demux/assignment/QC tables and consensus FASTA hashes.

Machine: macOS 26.6.1, ARM64, 16 logical CPUs, Python 3.12.10,
edlib 1.3.9.post1, PyYAML 6.0.3. Reporting extras were installed; their versions
are recorded with the measurements. Median of three runs:

| Requested execution | Total time | Input reads/s | Sampled peak process-tree RSS |
|---|---:|---:|---:|
| Serial, 1 worker | 3.252 s | 3,690 | 250.8 MiB |
| Thread, 4 workers | 3.373 s | 3,558 | 245.8 MiB |
| Process, 4 workers | 1.718 s | 6,984 | 271.0 MiB |

Four processes were **1.89× faster** than serial. Threads were about 4% slower.
All four compared scientific artifacts matched across **all nine runs**.
The implementation executes assignment serially for a requested thread backend;
this row measures that actual behavior, not four-thread assignment in isolation.

These are short repeated synthetic reads and a tiny reference panel, not a
production forecast. Timings include process startup and the installed reporting
environment; cache effects and system load matter. Memory is sampled every 20 ms
and sums RSS across parent and descendants, so shared pages can be counted more
than once and short peaks can be missed. It is not unique physical memory or a
whole-machine maximum. On restricted macOS environments, psutil process enumeration
can require broader execution permissions. Benchmark on the intended host.

The earlier 2026-08-12 lab entry measured real 20k-read data: demux approximately
591 reads/s serial versus 1,818 with eight processes; assignment approximately
660 versus 2,127 reads/s after accounting for startup. Those different workloads
support using processes but cannot be substituted for these fresh measurements.

### Configured parallelism versus this initial comparison

At the time of the sweep, production run profiles requested **eight process workers**, with
one native thread per worker and 1,000-read demultiplexing chunks. Profiles
using the `ont-r10-amplicon` preset inherit its process backend. `jobs: 0`
requests all detected CPUs (16 on this machine), subject to the resource planner.
Neither setting parallelizes every stage: consensus groups, QC and reporting
remain mostly serial.

The initial four-process comparison above did **not** test the production worker
count, maximum parallelism, or an optimum. Its 12,000 reads also supply only 12
demultiplexing chunks, too few to occupy 16 workers simultaneously. It established
output equivalence and a limited speedup; it should not have been the only fresh
parallel-performance evidence in the assessment.

Earlier real-data measurements favored eight workers: assignment was 2,127
reads/s with eight versus 2,021 with 16, and demux was 1,818 with eight versus
1,784 with 12. This explains the profiles' choice, but is not proof that eight
is best for every workload. Higher worker counts trade startup, index memory,
serialization and serial-stage limits against more concurrent computation.

### Expanded worker sweep — 120,000 synthetic reads

After the user challenged the four-process scope, repeated the comparison with
ten times as many reads (120 demux chunks), three trials per mode, and a seeded
shuffle of trial order. Explicitly included production's eight workers, all 16
CPUs, and `jobs: 0`. Source hashes, raw timings, memory samples, failure details
and output digests are in [the expanded evidence](benchmarks/2026-09-07-parallel-sweep.json).

```bash
python scripts/benchmark_portable.py \
  --output runs/parallel-sweep-20260907-complete --reads 120000 --repeats 3 \
  --process-jobs 1 2 4 8 12 16 0
```

Median full-run results on the same machine:

| Backend / requested workers | Resolved worker budget | Time | Reads/s | Peak tree RSS |
|---|---:|---:|---:|---:|
| Serial / 1 | 1 | 27.209 s | 4,410 | 292.7 MiB |
| Thread / 4 | 4 (assignment serial) | 29.022 s | 4,135 | 290.0 MiB |
| Process / 1 | 1 | Failed all 3 trials | — | — |
| Process / 2 | 2 | 16.227 s | 7,395 | 324.0 MiB |
| Process / 4 | 4 | 10.153 s | 11,819 | 443.0 MiB |
| Process / 8 | 8 | 7.719 s | 15,547 | 694.7 MiB |
| Process / 12 | 12 | 7.024 s | 17,085 | 927.5 MiB |
| Process / 16 | 16 | 6.937 s | 17,299 | 1,104.8 MiB |
| Process / 0 (detect) | 16 | 6.981 s | 17,189 | 1,120.7 MiB |

Sixteen processes were 3.92× faster than serial, 32% shorter in elapsed time
than four, and **10% shorter than eight**, with **59% more summed peak RSS than
eight**. Twelve was only 1.3% longer than 16; their individual timing ranges
overlap, so three trials do not establish a robust optimum between those counts.
Explicit 16 and automatic detection behaved similarly. This workload benefits
from using more than four workers, but remains a tiny-panel synthetic workload.
The result does not supersede representative production measurements.

All four compared scientific artifacts matched across **24 successful trials**.
All three one-process trials failed at assignment with `KeyError: 'config'`:
`_ordered_map` takes its serial shortcut without invoking the initializer, but
the caller still chooses `_assign_batch_worker`, which requires initialized
worker globals. This also affects a process configuration clamped to one CPU.
The benchmark now records failures, preserves incremental `trials.json`, and
continues the remaining modes. Its exit status is intentionally nonzero for this
sweep; output equivalence applies to successful runs, not the failed mode.

The sweep itself changed no production settings; the subsequent user request
changed the default and shipped profiles to automatic detection. MCP's default `--max-jobs 1`
limits concurrent
pipeline runs, not worker processes inside a run; it passes the YAML's parallel
settings through to the CLI. Maximum process parallelism is selected with
`parallel.backend: process`, `parallel.jobs: 0` and `threads_per_job: 1`, subject
to the existing planner limitations. The one-CPU defect above is now fixed.

### Follow-up after the four priority fixes

All **439 tests passed** on macOS/Python 3.12.10; the rebuilt wheel was installed
outside the checkout and completed the real stdio MCP example workflow. Its
installed source identity matches the checkout. Linux cgroup/affinity behavior
is covered by simulated v1/v2 fixtures; this session did not execute remote
Linux/Windows/Colab runners.

Three new 120,000-read trials with automatic **16-worker** process execution took
7.365, 7.429 and 7.584 seconds (median **7.429 s**). All nine checked artifacts
match the earlier serial run exactly after decompression, including contributor
IDs and consensus summary counts. The shipped seven-read example also matches
its saved pre-change baseline. Raw results, hashes and reproduction commands are
in [priority-fix verification](benchmarks/2026-09-07-priority-fixes.json).

The new median is about 6.4% longer than the earlier 6.981-second auto-worker
median. These are historical, unpaired timings with different monitoring
methods, so they indicate the scale of the overhead rather than isolate its
cause. RSS was not collected in this follow-up because this sandbox restricts
process-tree enumeration. A separate regression streams 8 MiB of synthetic
payload across 2,000 disk-backed groups with less than 2 MiB of traced Python
allocations; that excludes SQLite's native cache and is not an RSS measurement.

### Where time and memory go

- **Input I/O:** `validate_inputs` hashes the whole FASTQ and normally scans it
  again; demultiplexing reads it a third time. A separate `validate --quick`
  followed by `run` adds another complete hash pass. Quick means no record scan,
  not a cheap metadata check. Stage hashes add intermediate-file I/O too.
- **Parallel work:** edlib calls coexist with substantial Python dictionary,
  k-mer and serialization work. Explicit `backend: process` accelerates demux and
  assignment. `auto` does not automatically choose the process backend.
- **Scheduling:** `_ordered_map` keeps only one pending batch per worker and
  waits in input order. It bounds queued batches but can idle workers behind a
  slow early batch. Increase buffering only after measuring representative data.
- **Consensus/chimera analysis:** the follow-up now partitions eligible reads and
  routed wells on local SQLite storage and processes one group at a time. Claimed
  IDs and output accumulators also stay on disk. This adds disk/serialization work
  while removing retention across the entire dataset. The configurable per-group
  estimate guard fails explicitly when necessary; it is not a total RSS cap.
  Consensus groups still run serially.
- **Native threads:** both workers and threads are now clamped to the effective
  allocation, and the bounded thread count reaches `build_reference_consensus`.
- **Reports/storage:** gzip JSONL is portable and inspectable but duplicates
  sequences in intermediates. Per-clone output trees create many small files,
  expensive on Drive/network mounts. A cold or unwritable matplotlib font cache
  can add substantial startup time; MCP uses a writable workspace cache unless
  `MPLCONFIGDIR` is already supplied.

`scripts/benchmark_assignment.py` does not currently resolve/apply full-length
flanks and its serial path does not reproduce all pipeline gate/layout wiring.
Use whole-pipeline timing and per-read equivalence for full-length profiles until
that script shares the same preparation code. Matching status totals alone is
weaker evidence than comparing per-read decisions and final artifacts.

## Portability

| Target | Current position | Remaining qualification |
|---|---|---|
| macOS ARM64 | Core and MCP exercised here, including spawn workers | Production scale and optional native tool behavior still need workload-specific checks |
| macOS Intel | Standard Python/path/subprocess design | Not run on Intel in this assessment; verify dependency wheels |
| Linux x86-64 | Core CI configured for Python 3.10/3.13; spawn-safe workers, headless plots | Run updated CI; account for cgroups/affinity and memory limits |
| Native Windows x86-64 | Core CI configured for Python 3.12; `pathlib`, spawn, shell-free commands | Run MCP cancellation tests there; keep paths short; use portable consensus |
| WSL2 | Practical environment for Unix bioinformatics tools | Install/run Python and server inside WSL; prefer its Linux filesystem for active work |
| Linux/Windows ARM, Alpine | Conditional | edlib and optional compiled dependencies need matching wheels or a compiler/toolchain; not validated |
| Google Colab | Thin package-driven notebook, local staging guidance, optional MCP client cell | Hosted runtime not tested; transient disk, varying resources, supported Python version and Drive I/O matter |

The core is **not pure Python**: edlib is a compiled extension distributed with
platform wheels. Standard platforms can install easily, but “pip-installable”
does not guarantee a wheel for every architecture/interpreter. See
[edlib distributions](https://pypi.org/project/edlib/). Python 3.14 is explicitly
excluded by the project's `requires-python`; the MCP extra does not change that.

Native MAFFT/SPOA remain opt-in. The portable workflow does not need conda,
Bioconda, apt or those executables. No GPU-specific execution exists, so a Colab
GPU allocation does not directly accelerate this pipeline.

For Colab, copy FASTQ to `/content`, write intermediates locally, and archive
completed evidence to Drive. Do not promise recovery after VM recycling without
a durable backup; the current notebook's small result export does not preserve
all intermediates required for rerun. Archive the whole report stage if you need
its linked figures, not just `report.html`. A kernel-side MCP client can use stdio;
a remote desktop/cloud AI cannot reach the VM's loopback service directly.
Hosted resources and session lifetimes are variable, not guaranteed.
See the [official Colab FAQ](https://research.google.com/colaboratory/faq.html).

There is no committed environment lock or container image defining a canonical
scientific runtime. `environment.yml` uses an editable relative pip install,
so create it from the checkout. Git/Colab installs tracking a branch can change.
Pin a source revision plus platform-specific package/native-tool versions for
study-level reproducibility; test wheel installation outside the checkout.

## Prioritized findings

| Priority | Finding and evidence | Recommended next change |
|---|---|---|
| Fixed | Subsample could erase its own input or leave partial output | Reject existing destinations/aliases; validate in a temporary sibling and atomically create without clobbering; race and malformed-input regressions |
| Fixed | Zero consensus cap crashed and selected no chimera members | Zero includes all eligible reads; portable API, pipeline contributors and chimera output regressions |
| Fixed | Oversized native-thread requests exceeded the CPU budget | Clamp workers and native threads; pass bounded threads into consensus/MAFFT |
| Fixed | One process worker skipped initialization | Prepared local callable; fresh-interpreter regressions for requested 0, 1 and 8 workers on one CPU |
| Fixed | Host CPU count ignored process/container allocation | Effective process CPU count, affinity and visible cgroup v1/v2 quotas/cpusets, including ancestors |
| Fixed | Unchanged version number allowed stale code reuse | Installed source/preset content enters fingerprints/run metadata; legacy or different identities refuse reuse before copying |
| Fixed | All-well/all-group accumulation grew with dataset size | SQLite partitions groups, claimed IDs and output accumulators; per-group memory estimate guard defaults to 512 MiB and honors cgroup memory limits |
| Medium | Core CLI/API run IDs are concatenated directly into paths, and rerun can remove its new destination on failure | Share Windows-safe run-ID and destination validation with the CLI/API; MCP already validates its IDs |
| Medium | Assignment benchmark can diverge from full-length preparation and gates | Reuse a public prepared-stage context rather than rebuilding a partial context |
| Remaining | Stages after assignment are largely serial; reference indexes and worker batches still consume RAM | Profile representative runs before further parallelization; group admission guards are not total RSS limits |
| Medium | No lock/container; tests run primarily against source installations | Add reproducible release environments and observe the new installed-wheel CI job |
| Low | Historical architecture layout/vocabulary and the notebook's pure-Python claim were stale | Current layout now called out; notebook dependency wording corrected |

The fingerprint issue is also documented in the 2026-08-27 lab history: recording
a Git revision in `run.json` does not make it part of a stage's reuse decision.
Rerun's checksum verification confirms artifact integrity, not compatibility
between two different implementations that both call themselves `0.3.0`.

The subsample data-loss finding was reproduced during the follow-up priority
review using an isolated temporary FASTQ, not experimental data. It raises file
protection above structural refactoring in the cleanup order. The implementation
also leaves partially written output after a parsing failure. This finding is
preserved as historical evidence; file protection is now implemented and regression-tested.

Historical reproduction of the uncapped-consensus defect, using only shipped data
(run from the checkout; outputs are confined to a temporary directory):

```python
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
from nanopore3.config import load_config
from nanopore3.pipeline import run_pipeline

config = load_config("configs/example.yaml")
config = replace(config, consensus=replace(config.consensus, maximum_reads=0))
with TemporaryDirectory() as directory:
    run_pipeline(config, output_root=Path(directory), run_id="zero")
# Before the fix: IndexError in stage 04_consensus.
# Current behavior: completes using every eligible read (subject to the memory guard).
```

Python documents the distinction between system CPUs and CPUs usable by the
process; affinity-aware counts still require separate quota handling for some
container setups. See [Python CPU-count documentation](https://docs.python.org/3/library/os.html#os.process_cpu_count).

## Changes delivered with this assessment

- Optional MCP extra and executable; 12 tools, two resources and one prompt.
- Background CLI jobs with bounded concurrency, durable status/logs, cancellation,
  restart-aware reporting, bounded previews and workspace path validation.
- Local stdio and loopback Streamable HTTP; platform/client/Colab instructions and
  a runnable SDK smoke client. No public deployment or AI-account configuration.
- Portable synthetic benchmark and versioned measurements; no private input data.
- MCP tests added to the existing Linux/macOS/Windows matrix after core-only tests;
  a separate Linux wheel smoke job includes reporting and MCP extras.
- Build-system minimum raised from setuptools 69 to 77.0.3 because the existing
  SPDX license string and `license-files` require the newer metadata support.
  See [setuptools licensing documentation](https://setuptools.pypa.io/en/stable/userguide/license_migration.html).
- README/architecture links and notebook corrections; no analysis algorithms,
  experimental configuration, thresholds or existing outputs changed.

Remaining follow-up: share core CLI/API run-ID validation, lock release
environments, and profile representative large inputs before refactoring stages
or adding more parallelism. The four requested fixes are covered by
`tests/test_priority_fixes.py` plus existing pipeline/rerun/backend tests.
