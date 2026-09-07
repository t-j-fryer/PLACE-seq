# Nanopore3

Nanopore3 is a clean-room successor to the exploratory Nanopore2 notebooks. Its goal is a
portable, reproducible, and inspectable pipeline for demultiplexing Nanopore amplicons,
assigning reads to references, constructing consensuses, and reporting quality-control
evidence.

> **Status: v0.3 pre-alpha.** The current release implements a portable end-to-end workflow,
> explicit multi-library plate routing, and the contracts needed for experimental validation. It is not yet a
> production-validated replacement for Nanopore2. Keep Nanopore2 and its results unchanged
> while outputs are compared against synthetic controls and held-out experimental runs.

**September 2026 update:** maximum allocated parallelism, safe subsampling,
source-aware reuse, disk-backed grouping and a documented MCP server. Read the
[update and migration guide](docs/updates/2026-09-07-portability-and-mcp.md)
before resuming an older run.

## Design promises

- The portable core installs with pip on macOS, native Windows, Linux, and Google Colab.
- Inputs are never modified. A run writes to a new directory and records its configuration.
- Classification has explicit rejection states; weak evidence is not silently called assigned.
- Parallel execution is bounded, deterministic, and avoids nested CPU oversubscription.
- Optional native tools are detected during preflight and recorded in provenance.
- Notebooks are clients of the package, not the implementation or source of hidden state.

**Try it without installing anything:**
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/t-j-fryer/nanopore3/blob/main/notebooks/Nanopore3_Colab.ipynb)
— runs the synthetic example after Setup, then walks through your own data from Drive.
For private repository access, use an authenticated checkout or upload a wheel.

**New here? Start with [Worked examples](docs/workflows.md)** — three complete scenarios,
from what is on the bench to what you read afterwards.

**Using an AI assistant?** The optional [MCP server](docs/mcp.md) exposes project
inspection, configuration creation, safe FASTQ subsampling, validation, background runs, resume/rerun,
job logs, cancellation and result previews through stdio or Streamable HTTP:

```bash
python -m pip install -e ".[mcp]"
python -m nanopore3.mcp_server --workspace /absolute/path/to/project
```

The guide includes Windows setup, client configuration, an executable smoke
client, and Colab usage. See the [repository assessment](docs/repository-assessment.md)
for measured performance, portability evidence and known core limitations.

See [Architecture](docs/architecture.md) for the stage model and reproducibility contract,
[Configuration](docs/configuration.md) for the option reference,
[References](docs/references.md) for the four ways to describe your constructs, and
[Pooling layout](docs/pooling-layout.md) if several culture plates were combined into one
colony-PCR plate.

## Quick start

Nanopore3 requires Python 3.10–3.13. Python 3.12 is the reference development version.

Once installed, create and run a self-contained example:

```bash
nanopore3 init nanopore3-example
nanopore3 validate --config nanopore3-example/configs/example.yaml
nanopore3 run --config nanopore3-example/configs/example.yaml
```

Each run gets a new directory containing checksummed stage manifests, per-read decisions,
consensus contributor IDs, tri-state QC and an HTML report. Use `--run-id NAME --resume` only to
resume that exact configuration; changed or corrupted artifacts are rejected.

### macOS and Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
nanopore3 --help
```

### Windows PowerShell

```powershell
py -3.12 -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
nanopore3 --help
```

Native Windows supports the portable edlib-based workflow. WSL2 is recommended when a run
requires Unix bioinformatics binaries such as MAFFT, SPOA, or minimap2.

### Micromamba or Conda

The supplied environment is deliberately portable and does not require Bioconda:

```bash
micromamba create -f environment.yml
micromamba activate nanopore3
nanopore3 --help
```

### Jupyter and Google Colab

The [notebook](notebooks/Nanopore3_Colab.ipynb) works in local Jupyter on
Linux/macOS/Windows and in Colab. Clone or upload this repository, then install
into the active Colab kernel:

```python
%pip install -e "/content/Nanopore3[notebook,report,mcp]"
```

Notebook Setup prefers a nearby checkout, or accepts an explicit wheel/pinned
package specification, and installs reporting/MCP together.
The notebook uses checked CLI subprocesses through the kernel's interpreter,
shows effective CPU/memory budgets, and creates new paths for inputs and runs.
Stage active data/output on local disk; archive to Drive or other durable storage.
Keep a full run archive for reuse, and reinstall the same source revision/wheel.
Summary exports alone cannot resume a run. The portable backend needs no apt or conda.

## Optional features

Install only what is needed:

```bash
python -m pip install -e ".[report]"          # tables, statistics, and plots
python -m pip install -e ".[notebook]"        # JupyterLab and kernel support
python -m pip install -e ".[report,notebook]"
python -m pip install -e ".[dev]"             # tests, coverage, and linting
```

MAFFT and SPOA provide the opt-in `mafft_spoa` consensus backend used by the legacy
notebook; `portable` is the default reference-guided edlib pileup. Minimap2 and mappy are
reserved optional accelerators. These tools are intentionally not installed by the core package
because availability differs across operating systems. Nanopore3 locates requested executables,
reports their versions, and fails preflight when a selected backend is unavailable. The actual
consensus backend is recorded with each result.

For scientifically locked production runs, use a platform-specific environment lock or a
versioned Linux container in addition to the portable package metadata.

## Commands

```bash
nanopore3 doctor --json
nanopore3 validate --config configs/example.yaml
nanopore3 run --config configs/example.yaml --output runs
```

The same workflow is available to Python and notebooks through `load_config()` and
`run_pipeline()`. See [Configuration](docs/configuration.md) and the clean quick-start notebook.

Barcode panels live in reviewable CSV registries and run profiles select one primer family
explicitly. Named reference libraries are routed by plate barcode rather than pooled. See
[Barcode registry](docs/barcode_registry.md) and the
[20260506 pilot report](docs/20260506_lab_biotin_pilot.md). The recovered
thresholds, decision rules and parallel benchmarks are documented in the
[demultiplex optimisation audit](docs/demultiplex_optimisation.md) and the
[assignment optimisation note](docs/assignment_optimisation.md).

## Lab notebook

[`LAB_NOTEBOOK.md`](LAB_NOTEBOOK.md) is the dated, append-only record of what changed in this
repository, why, what was learned, and what to do next. **Read the most recent entry before
starting work, and add an entry when you finish.** It is written so that a person or an agent with
no prior context can pick up the work. Entries separate scientific changes, which can alter
results, from performance changes, which must not.

## Performance

By default, `parallel.backend: process`, `parallel.jobs: 0`, and
`parallel.threads_per_job: 1` use all allocated CPUs for demultiplexing and
assignment. Shipped profiles use these settings too. Set a positive `jobs`
value to cap workers, or `backend: serial` to disable process parallelism.
Linux affinity/cgroup limits are respected. Consensus/chimera groups spill to
disk and have a configurable memory estimate guard; `maximum_reads: 0` includes
all eligible reads. See the [resource and reuse guide](docs/configuration.md).
Resume/rerun require unchanged installed source; older runs without source identity
remain readable and need a fresh run for recomputation.

Reads are matched to references with a k-mer prefilter before any alignment: a specificity-weighted
index proposes a few candidates and only those are aligned, which costs about four alignments per
read instead of one per reference. Both demultiplexing and assignment are bound by Python-level work
rather than by the GIL-releasing edlib calls, so **process workers scale this pipeline and threads do
not** — on a 16-CPU Mac, four threads demultiplexed *slower* than serial, while eight processes were
roughly three times faster than serial. Benchmark your own machine with `scripts/benchmark_demux.py`
and `scripts/benchmark_assignment.py` before changing `parallel` in a run profile, and confirm the
call counts are unchanged.

## Repository layout

```text
LAB_NOTEBOOK.md    dated record of changes, rationale, lessons, and next steps
src/nanopore3/     installable library and CLI
configs/           versioned run and assay configuration examples
scripts/           benchmarking and diagnostic utilities
fixtures/          small synthetic data suitable for version control
tests/             unit, integration, and golden tests
docs/              architecture and operating guidance
notebooks/          package-driven local/Colab examples (no workflow implementation)
```

Large FASTQ data and generated runs are ignored by Git. Reference manifests and checksums should
be versioned; raw sequencing data should be managed by the laboratory's durable data system.

## Safety

Do not point experimental development at the only copy of a dataset. Nanopore3 treats inputs as
read-only and should write intermediate files atomically. It must not infer completion merely
because some output files exist, and it must not delete or move prior evidence during analysis.
Nanopore3 includes opt-in positional-signature chimera analysis; its reported origin
classifications remain dependent on configured references and experimental validation. Until
the experimental validation suite is complete, compare every result with controls and retain the
original Nanopore2 analysis.
