# Nanopore3

Analyse Nanopore amplicon reads from FASTQ to a browsable report. Nanopore3
separates reads by plate and well barcode, matches them to reference sequences,
builds consensus sequences, and reports quality control (QC) for each result.
Use it from the command line, a Jupyter/Colab notebook, or an AI assistant with
Model Context Protocol (MCP) support. All three use the same pipeline.

> **v0.3.0 · pre-alpha.** The workflow runs end to end, but experimental
> validation is still in progress. Compare results with controls and retain
> existing Nanopore2 analyses before adopting it as a replacement.

## Choose how to run

| I want to… | Start here |
|---|---|
| Run on my computer or a Linux server | [Install and run the example](#install-and-run-the-example) |
| Work interactively in Jupyter or Google Colab | [Notebook setup](#jupyter-and-google-colab) |
| Let an AI assistant run and inspect analyses | [Connect an AI assistant](#connect-an-ai-assistant) |
| Configure my own sequencing experiment | [Use your own data](#use-your-own-data) |

**Updating an existing installation?** Read the
[update and migration guide](docs/updates/2026-09-07-portability-and-mcp.md)
before resuming older runs. Reuse requires matching installed source code.

## Install and run the example

You need **Python 3.10–3.13**; use Python 3.12 if you are setting up a new
environment. The portable workflow runs on macOS, Linux and native Windows
without additional bioinformatics executables.

### 1. Get the repository

With Git installed and access to this repository:

```sh
git clone https://github.com/t-j-fryer/nanopore3.git
cd nanopore3
```

Alternatively, use GitHub's **Code → Download ZIP**, extract it, and open a
terminal in the extracted folder containing `pyproject.toml`. If you already
have a checkout, open that folder. This update is distributed through the
repository or a wheel built from it; it is not a PyPI release.

### 2. Install in a virtual environment

**macOS / Linux**

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[report]"
nanopore3 doctor
```

Use `python3` instead of `python3.12` if it points to a supported Python version.

**Windows PowerShell**

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[report]"
nanopore3 doctor
```

If PowerShell blocks activation, use `.\.venv\Scripts\python.exe` instead of
`python` in the install commands, and `.\.venv\Scripts\python.exe -m nanopore3`
instead of `nanopore3` in subsequent commands. No execution-policy change is needed.

The `report` extra adds plots and their dependencies. For a smaller installation,
use `python -m pip install -e .`; consensus sequences, core tables and the HTML
report still work. Missing optional executables in `doctor` are expected when
using the default `portable` consensus backend.

<details>
<summary>Prefer Micromamba or Conda?</summary>

From the repository folder:

```sh
micromamba create -f environment.yml
micromamba activate nanopore3
nanopore3 doctor
```

For Conda, replace `micromamba` with `conda env` in the create command and use
`conda activate nanopore3`. The environment includes reporting and notebook
support; it does not require Bioconda.

</details>

### 3. Run the bundled synthetic example

Run these commands from the same terminal. No sequencing data is needed:

```sh
nanopore3 init nanopore3-example
nanopore3 validate --config nanopore3-example/configs/example.yaml
nanopore3 run --config nanopore3-example/configs/example.yaml --run-id first-run
```

The CLI prints stage progress and finishes with `Completed run:` and the output
paths. Results are under `nanopore3-example/runs/first-run/`. To run the example
again, choose a new `--run-id`, such as `second-run`.

### 4. Open the results

Paths below are relative to the completed run directory:

| File or folder | What to look at |
|---|---|
| `stages/06_report/report.html` | Open in a browser for the run overview |
| `consensus_by_plate/index.csv` | Open in a spreadsheet for one row per clone, grades and sequence accuracy |
| `consensus_by_plate/` | Individual consensus FASTA files organised by plate and well |
| `stages/05_qc/qc.csv.gz` | Detailed QC calls in a compressed CSV |
| `run.json` and stage manifests | Recorded configuration, source identity, checksums and execution details |

Clone files and their index are written by default (`consensus_tree: true`).
A completed command means the analysis finished; inspect QC to determine which
results passed, failed or remain uncertain. See
[Reading the results](docs/workflows.md#reading-the-results) for column meanings.

## Use your own data

Start from a generated example, then edit its YAML configuration to describe
**your experiment**. The bundled synthetic motifs, barcodes and thresholds are
for the demo.

You will need:

- **Reads:** one or more FASTQ or FASTQ.gz files.
- **References:** FASTA sequences for your designed inserts or assembled constructs.
- **Assay details:** boundary motifs, plate/well barcode sequences or CSV registries,
  and the reference library associated with each plate barcode.
- **Pooling layout, if applicable:** the mapping from pooled PCR plates back to
  source culture plates.

The [four worked examples](docs/workflows.md) cover inserts versus assembled
references, with one culture plate per barcode or pooled plates. Use the
[configuration guide](docs/configuration.md) for individual settings.
**Relative file paths resolve beside the YAML file**, not beside your terminal.

After saving your configuration as `run.yaml`:

```sh
nanopore3 validate --config run.yaml
nanopore3 run --config run.yaml
```

Input files are read-only. New runs receive their own output directory under
`output_root`; `--output` overrides that location. Keep original reads and
complete run directories in durable storage.

## Jupyter and Google Colab

The same [notebook](notebooks/Nanopore3_Colab.ipynb) works in local Jupyter on
macOS, Linux and Windows, and in Colab.

**Local Jupyter:** after the installation above, run from the repository folder:

```sh
python -m pip install -e ".[report,notebook,mcp]"
jupyter lab notebooks/Nanopore3_Colab.ipynb
```

**Google Colab:** open the notebook below, or upload the `.ipynb` file to Colab.
No local Python installation is required.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/t-j-fryer/nanopore3/blob/main/notebooks/Nanopore3_Colab.ipynb)

**Private repository access:** the badge and the notebook's default Git install
require access to the repository. If either fails, upload the notebook and a
wheel built from an accessible checkout. Build the wheel locally with
`python -m pip wheel . --no-deps --wheel-dir dist`, upload the resulting `.whl`
to Colab, and set `PACKAGE_SPEC` in Setup to its path with `[report,mcp]` appended
(for example, `/content/nanopore3-0.3.0-py3-none-any.whl[report,mcp]`).

Run **Setup**, then choose **A — Synthetic example** or **B — Your own data**.
The notebook includes data staging, resource settings, results and archiving.
In Colab, keep active inputs and outputs on VM-local disk, then archive to Drive.
The VM's files are temporary; retain a **complete run archive** for resume/rerun.
Live Colab execution has not yet been verified; local notebook workflows and a
simulated Drive mount are tested.

## Connect an AI assistant

The optional MCP server lets a compatible AI host inspect a project, create
configs, subsample reads, validate inputs, run analyses, monitor or cancel jobs,
and read results. The host must support **stdio or Streamable HTTP MCP**.

From the installed checkout:

```sh
python -m pip install -e ".[report,mcp]"
python -m nanopore3.mcp_server --help
```

Follow the [MCP setup guide](docs/mcp.md#connect-a-local-ai-using-stdio) to configure
your host with the environment's **absolute Python path** and an **existing
workspace directory**. It includes copyable JSON, Windows paths, all 13 tools,
HTTP setup and a smoke client. For stdio, the AI host launches the server.

Suggested first request once connected:

> Create the Nanopore3 synthetic example, wait for it to finish, validate it,
> then run it. Monitor each job to completion and show me the report and QC results.

For Colab, the MCP client must run inside the notebook VM; a desktop client
cannot directly reach that VM's loopback server.

## Performance and resource settings

**Maximum allocated CPU parallelism is the default** for demultiplexing and
assignment:

```yaml
parallel:
  backend: process
  jobs: 0
  threads_per_job: 1
```

`jobs: 0` uses the detected CPU allocation, respecting exposed affinity and Linux
cgroup limits. Set a positive `jobs` value to use fewer workers, or
`backend: serial` for serial execution. Use `nanopore3 doctor --json` for runtime
information and `nanopore3 validate --config run.yaml --quick` for the configured
resource plan; even quick validation computes input checksums.

More workers also need more RAM. Grouped reads use temporary disk storage with
an estimated per-group memory guard; this is not a total process memory limit.
`consensus.maximum_reads: 0` includes all eligible reads in each consensus;
positive values cap depth. Maximum CPU use does **not** change that read cap.
See [resource settings](docs/configuration.md#running-on-a-machine-you-did-not-configure-for)
and the [performance assessment](docs/repository-assessment.md) for controls,
measured benchmarks and their limitations.

## Common tasks and troubleshooting

| Task or problem | Command or next step |
|---|---|
| Explore commands | `nanopore3 --help` or `nanopore3 run --help` |
| Check installation and CPUs | `nanopore3 doctor` |
| Try the first 100,000 reads | `nanopore3 subsample --input reads.fastq.gz --output subset.fastq.gz --reads 100000` — use a new output path on local disk; this selects a prefix, not a random sample |
| Inspect a pooling layout | `nanopore3 layout --config run.yaml` |
| Resume an interrupted run | `nanopore3 run --config run.yaml --run-id NAME --resume` — use the original output root, configuration and installed source |
| Recompute selected stages | See the [rerun guide](docs/architecture.md#reruns-recomputing-the-analysis-without-the-fastq) |
| `nanopore3` command not found | Activate the environment, or use its Python executable with `-m nanopore3` |
| Missing files or environment variables | Check paths relative to the YAML and set every `${VARIABLE}` it uses |
| Source identity missing or changed | Old results remain readable; start a new run from original inputs. Restart Python/MCP after changing installed code |
| Group memory estimate exceeded | Review the [memory controls](docs/configuration.md#running-on-a-machine-you-did-not-configure-for) and available RAM/disk before increasing the limit |

The optional `mafft_spoa` consensus backend requires separately installed MAFFT
and SPOA. Use Linux or WSL2 on Windows for that toolchain; the default portable
backend needs neither. See the [backend policy](docs/architecture.md#portability-and-backend-policy).

## Documentation and development

| Guide | Contents |
|---|---|
| [Worked examples](docs/workflows.md) | Choose an assay configuration and interpret results |
| [Configuration](docs/configuration.md) | Paths, presets, resources, QC and reuse |
| [References](docs/references.md) | Describe inserts, constant regions and assembled constructs |
| [Barcode registry](docs/barcode_registry.md) · [Pooling layout](docs/pooling-layout.md) | Barcode panels and source-plate routing |
| [MCP server](docs/mcp.md) | AI host setup, tools, jobs and troubleshooting |
| [Repository assessment](docs/repository-assessment.md) | Structure, performance, platform evidence and limitations |
| [Architecture](docs/architecture.md) | Pipeline stages and reproducibility contracts |
| [Update and migration guide](docs/updates/2026-09-07-portability-and-mcp.md) | September 2026 changes and compatibility |
| [Lab notebook](LAB_NOTEBOOK.md) | Dated changes, validation evidence and next steps |

Implementation lives in `src/nanopore3/`; examples and profiles in `configs/`;
synthetic data in `fixtures/`; checks in `tests/`; utilities in `scripts/`;
and the notebook in `notebooks/`. Large reads and generated runs are Git-ignored.

Before contributing, read the latest [lab notebook](LAB_NOTEBOOK.md) entry and
add an entry when finished. Install development tools with
`python -m pip install -e ".[dev,mcp]"`, then run `python -m pytest -q` and
`python -m ruff check src scripts tests`. Scientific tuning history is recorded
in the [demultiplex audit](docs/demultiplex_optimisation.md),
[assignment note](docs/assignment_optimisation.md) and
[pilot report](docs/20260506_lab_biotin_pilot.md).

Licensed under the [MIT License](LICENSE).
