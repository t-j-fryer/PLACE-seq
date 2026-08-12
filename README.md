# Nanopore3

Nanopore3 is a clean-room successor to the exploratory Nanopore2 notebooks. Its goal is a
portable, reproducible, and inspectable pipeline for demultiplexing Nanopore amplicons,
assigning reads to references, constructing consensuses, and reporting quality-control
evidence.

> **Status: v0.2 pre-alpha.** The current release implements a portable end-to-end workflow,
> explicit multi-library plate routing, and the contracts needed for experimental validation. It is not yet a
> production-validated replacement for Nanopore2. Keep Nanopore2 and its results unchanged
> while outputs are compared against synthetic controls and held-out experimental runs.

## Design promises

- The portable core installs with pip on macOS, native Windows, Linux, and Google Colab.
- Inputs are never modified. A run writes to a new directory and records its configuration.
- Classification has explicit rejection states; weak evidence is not silently called assigned.
- Parallel execution is bounded, deterministic, and avoids nested CPU oversubscription.
- Optional native tools are detected during preflight and recorded in provenance.
- Notebooks are clients of the package, not the implementation or source of hidden state.

See [Architecture](docs/architecture.md) for the stage model and reproducibility contract.

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

### Google Colab

Clone or upload this repository, then install it into the active Colab kernel:

```python
%pip install -e "/content/Nanopore3[notebook,report]"
```

Run the same Python API or CLI from the notebook. Keep datasets on mounted Drive or Colab
storage rather than committing them to the repository. The portable backend must work without
apt or conda; an optional setup cell may install native tools for an enhanced run.

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
[20260506 pilot report](docs/20260506_lab_biotin_pilot.md).

## Repository layout

```text
src/nanopore3/     installable library and CLI
configs/           versioned run and assay configuration examples
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
Nanopore3 v0.2 intentionally does not make final biological chimera calls. It reports uncertain
reads conservatively while breakpoint-aware classification is developed and validated. Until
the experimental validation suite is complete, compare every result with controls and retain the
original Nanopore2 analysis.
