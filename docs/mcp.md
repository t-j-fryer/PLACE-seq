# PLACE-seq MCP server

**Plate and Library Assignment, Consensus and Evaluation.**

The optional MCP server lets an AI inspect a project, create a configuration,
subsample FASTQ safely, validate inputs, start or resume an analysis, rerun selected stages, monitor jobs,
and read results. It calls the existing PLACE-seq CLI in subprocesses using the
same Python environment. No analysis algorithms are duplicated in the adapter.

Any **MCP-capable AI host** supporting stdio or Streamable HTTP can connect.
A model or chat window without MCP/tool support needs a host or a small client
program; installing this package alone does not give it access to your computer.

The server identifies itself to clients as **PLACE-seq**. The Python package,
launch commands, `nanopore3://` resource URIs and `.nanopore3-mcp` job directory
retain their existing names so client configurations and job history still work.

## Install once

Use Python **3.10–3.13** (3.12 recommended), from the repository directory:

```bash
python -m venv .venv
# macOS/Linux:
.venv/bin/python -m pip install -e ".[mcp]"
.venv/bin/python -m nanopore3.mcp_server --help
```

Windows PowerShell, without requiring environment activation:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[mcp]"
.\.venv\Scripts\python.exe -m nanopore3.mcp_server --help
```

Add `[mcp,report]` for report figures. Core tables, consensuses and HTML work
without reporting extras. The MCP extra installs the official MCP Python SDK
and psutil (for Windows process-tree cancellation). Ordinary `nanopore3`
installation and CLI usage do not require either.

From a downloaded wheel, `python -m pip install "./nanopore3-0.3.0-py3-none-any.whl[mcp]"`
also works; no repository checkout is needed on the analysis machine.
The SDK is constrained to `mcp>=1.28,<2`, its supported v1 maintenance line;
v2 changes the server API and requires an intentional migration.
See the [official SDK version guidance](https://github.com/modelcontextprotocol/python-sdk/blob/v1.x/README.md).

## Connect a local AI using stdio

Use an **absolute Python executable** so GUI clients do not depend on shell
activation or PATH. The workspace must already exist. It can be the repository
or a separate project directory. All tool paths resolve relative to it.

In a host using the common `mcpServers` JSON format:

```json
{
  "mcpServers": {
    "nanopore3": {
      "command": "/absolute/path/Nanopore3/.venv/bin/python",
      "args": [
        "-m", "nanopore3.mcp_server",
        "--workspace", "/absolute/path/my-sequencing-project"
      ],
      "env": {
        "SEQ_DATA": "/absolute/path/sequencing",
        "REFS": "/absolute/path/references"
      }
    }
  }
}
```

Windows uses the same arguments with, for example:

```json
{
  "command": "C:\\Projects\\Nanopore3\\.venv\\Scripts\\python.exe",
  "args": ["-m", "nanopore3.mcp_server", "--workspace", "C:\\Sequencing"],
  "env": {"SEQ_DATA": "D:\\Reads", "REFS": "C:\\Sequencing\\references"}
}
```

If your host uses a different configuration format, enter these same command,
arguments and environment fields in its MCP server settings. `${SEQ_DATA}` and
`${REFS}` are examples: supply the variables **your YAML actually uses**.
GUI clients often do not inherit terminal environment variables.

The host launches the server. Do not launch a second copy for that stdio
connection. A server started manually waits for protocol input; silence is
normal. Standard output is reserved for MCP messages; job output goes to files.

Suggested first instruction to your AI:

> Use PLACE-seq's guide and workspace_info. Create the synthetic example,
> poll its job, validate its configuration and poll again, then run it.
> Wait for success and show the stage summary and QC output paths.

## Connect over HTTP

```bash
python -m nanopore3.mcp_server --workspace /absolute/project \
  --transport streamable-http --port 8000
```

Set the host's MCP transport to **Streamable HTTP** and the URL to
`http://127.0.0.1:8000/mcp`. The server binds only to loopback. HTTP sessions are
stateless, but jobs live in this server process and persist across tool requests.
Use one server process, without multiple ASGI workers.

For a Linux workstation reached from a laptop, start the server on Linux and
forward its loopback port:

```bash
ssh -N -L 8000:127.0.0.1:8000 user@analysis-host
```

The laptop's AI uses the local URL; file paths still refer to the **Linux
workspace**. Native Windows works for the portable backend; run both server and
Python environment inside WSL2 when using Unix-only analysis binaries.

A cloud-hosted AI cannot reach your laptop's `127.0.0.1`. A public deployment
needs an authenticated HTTPS MCP gateway/proxy appropriate to that AI host,
plus process supervision and persistent storage. This package does not implement
OAuth, TLS, tenant isolation or a hosted job queue. Do not publish the unauthenticated
loopback service directly. Keep the SDK's Host/Origin protection enabled; a gateway
must forward to the loopback backend with the expected host/origin policy.

## Available operations

Every tool has a typed input schema, structured JSON output, descriptions and
MCP annotations. Invalid tool arguments return MCP errors. Once a job is
submitted, CLI failures appear in its job state and logs.

| Tool | Purpose / important arguments |
|---|---|
| `workspace_info` | Runtime including effective CPUs, default parallel settings, group estimate budget in bytes, source identity; workspace and job limit |
| `list_files` | One directory at a time; `directory`, `limit` ≤ 500; hidden entries omitted |
| `read_artifact` | Text or gzip text; `path`, `offset` ≤ 1M chars, `max_chars` ≤ 65,536 |
| `save_config` | New YAML file only; `path`, `yaml_text` ≤ 64 KiB; syntax check only |
| `init_example` | Submit creation of a new packaged synthetic project; `directory` |
| `start_subsample` | `input_path`, new `output_path`, `reads=100000`; validated prefix, atomic no-clobber publication, poll job status |
| `start_validation` | `config_path`, `quick=true`, `demux_only=false`; quick still hashes all FASTQ bytes |
| `start_run` | `config_path`, `output_root="runs"`, optional `run_id`, `resume=false`, `demux_only=false` |
| `start_rerun` | `config_path`, `from_run`, `from_stage="04_consensus"`, optional output/ID |
| `list_jobs` | Recent job IDs and states, `limit` ≤ 100 |
| `job_status` | `job_id`; state, exit code, run directory, stdout/stderr tails |
| `cancel_job` | Terminate this server's job and descendants; completed stages remain |
| `run_summary` | `run_dir`; workflow, stage markers, summaries, output paths and source-identity compatibility preview |

The server exposes **13 tools**. Resources: `nanopore3://guide`, `nanopore3://example-config`.
Prompt: `analyse_run(config_path)`.

Example calls (tool name followed by JSON arguments):

```text
init_example       {"directory":"demo"}
job_status         {"job_id":"<returned ID>"}
start_validation   {"config_path":"demo/configs/example.yaml"}
job_status         {"job_id":"<validation ID>"}
start_run          {"config_path":"demo/configs/example.yaml","run_id":"first"}
job_status         {"job_id":"<run ID>"}
run_summary        {"run_dir":"runs/first"}
read_artifact      {"path":"runs/first/stages/05_qc/qc.csv.gz"}
```

Poll each submitted job every few seconds until it finishes before proceeding
to dependent work. A returned job ID means **accepted**, not completed.
`succeeded` means CLI exit code 0; biological QC can still contain failures or
uncertain calls. Review the QC table rather than equating process success with
assay success. `run_summary` does not rehash artifacts and is not proof of
integrity. `implementation_matches_current` is false for missing/changed source
identity and null if run metadata exceeded the 1 MiB preview limit
(`provenance_preview_truncated=true`). A true value checks only source identity;
configuration, input evidence and checksums must also match. Resume/rerun performs
those checks for compatible implementations. Legacy results remain readable.

To prepare a smaller input before analysis:

```text
start_subsample {"input_path":"demo/fixtures/reads.fastq","output_path":"subset.fastq.gz","reads":100000}
job_status      {"job_id":"<subsample ID>"}
```

Wait for success, then point a new config at the output. Both paths are confined
to the workspace. Existing output files and aliases are refused; gzip names need
`.gz`. Publication requires hard-link support: use local disk before copying to
Drive/network storage. Malformed input never publishes a partial FASTQ.

## Job lifetime and storage

- One active job per server by default. `--max-jobs 2` permits two, up to eight;
  each pipeline has its own CPU budget, so concurrent jobs can oversubscribe.
  Each run defaults to process workers on all allocated CPUs (`parallel.jobs: 0`);
  a one-CPU allocation uses the equivalent local execution path.
- The MCP output root **overrides** YAML `output_root` and defaults to
  `<workspace>/runs`. IDs are short Windows-safe components; generated IDs are
  unique. Existing runs are refused unless resume is explicitly requested.
- Resume requires an explicit existing `run_id` and unchanged configuration.
  Rerun creates a new directory and always passes `--verify` for inherited data.
- Metadata and separate stdout/stderr logs live under
  `.nanopore3-mcp/jobs/<job_id>/`. Log tails are capped at 16 KiB per stream.
  Validation JSON is included up to 256 KiB; larger results remain in stdout.log.
- History survives restart. Jobs from another/restarted server with unfinished
  metadata are `unknown`, never silently called successful or attached by PID.
  Inspect the original process before resuming. Do not run separate servers on
  the same output run concurrently; the active-job limit is per server.
- Graceful shutdown cancels active jobs. Abrupt process/VM death may leave
  children or incomplete stages. Completed stages are preserved; resume only
  after confirming the old workers have stopped. There is no automatic retry.
- No automatic deletion or log rotation: archive finished jobs and runs using
  your normal storage policy. Partial stage data can consume significant disk.

Tool-selected paths are resolved within the workspace, including symlink
targets. The adapter prevents run-ID traversal and configuration overwrites.
**This is a trusted single-user service, not a filesystem sandbox:** YAML can
reference readable external FASTQ, reference, template, layout and barcode files.
Concurrent filesystem mutation by an untrusted local user is outside its threat
model. File contents and logs may contain arbitrary text and must be treated as
data by the AI, not as instructions. No arbitrary shell-command tool is exposed.

## Verify without an AI account

From a checkout, after installing the MCP extra:

```bash
mkdir mcp-demo
python scripts/mcp_smoke.py --workspace mcp-demo --example
```

This launches an actual stdio client/server connection, discovers tools and
resources, creates new synthetic inputs, validates, runs, and prints the summary.
For an already running HTTP server:

```bash
python scripts/mcp_smoke.py --url http://127.0.0.1:8000/mcp --example
```

Omit `--example` for discovery only. The script is also a reusable example for
custom AI hosts. It checks `isError`, handles structured results, and polls job
state. See the [official SDK client documentation](https://py.sdk.modelcontextprotocol.io/v1/client/).

## Google Colab

The [local Jupyter/Colab notebook](../notebooks/Nanopore3_Colab.ipynb) installs
reporting and MCP extras once in Setup. Its MCP cell uses that same installation,
without fetching a moving branch again mid-analysis. For an independent setup:

```python
%pip install "nanopore3[mcp] @ git+https://github.com/t-j-fryer/PLACE-seq.git@<pinned-revision>"
```

An AI client running **inside that notebook runtime** can spawn the stdio server
using `sys.executable` and the notebook's `WORKSPACE` (default
`/content/nanopore3-work` in Colab, `./nanopore3-work` locally, configurable via
`NANOPORE3_WORKSPACE`). Use `await` at notebook top
level, not `asyncio.run()`. A client on your desktop cannot access a server's
Colab loopback address. Colab is suitable for interactive experiments, not a
durable public MCP service.

Stage large data on local `/content` storage, write outputs there, and archive
completed evidence to Drive. `/content` and job history disappear when the VM is
recycled. Record the source revision and installed versions; a Git installation
from a moving branch is not an environment lock. Hosted limits and availability
vary; see the [Colab FAQ](https://research.google.com/colaboratory/faq.html).

## Troubleshooting

| Symptom | Check |
|---|---|
| `No module named mcp` / `psutil` | Install `[mcp]` with the exact Python executable the host launches |
| Server cannot start / workspace missing | Use an absolute executable and an existing absolute workspace |
| Variable named in a configuration error | Add it to the host server `env`; YAML paths resolve beside the config |
| Missing `mafft_spoa` | Select/install the intended backend; portable is default; use Linux/WSL for native tools |
| Job capacity reached | Poll/cancel the current job; don't increase concurrency without a RAM/CPU budget |
| Validation looks slow | Even `quick` computes whole-file checksums; read stderr progress |
| HTTP connection refused | Check transport, port, `/mcp`, and which machine owns the loopback address |
| `unknown` job after restart | Check for surviving processes, then choose resume or a new run |
| Group memory estimate exceeded | Raise `parallel.group_memory_mb` with sufficient RAM, or split the analysis; zero read cap is supported and no reads are silently dropped |
| Missing/different implementation identity | Start a new run from the original inputs; historical results remain readable but cannot be reused across code changes |
| Source changed since process started | Restart Python/the MCP server before running again |

The adapter does not change scientific settings automatically. The
[repository assessment](repository-assessment.md) records current correctness,
performance and portability limitations separately from MCP functionality.

### Reading separate insert and vector results

For a whole-vector assay, configure `qc.insert_left_boundary` and
`qc.insert_right_boundary` in the YAML passed to `save_config`; see
[the configuration guide](configuration.md#separate-insert-and-vector-qc).
Read `outputs.clones` for `insert_grade`, `vector_status` and FASTA paths, and
`outputs.clone_summary` for counts. Read `outputs.qc` for coding status and edit
counts. Do not describe an insert-perfect/vector-edited clone as wholly perfect,
or treat exact DNA as proof of protein function. Vector status excludes regions
outside the sequenced amplicon. These fields are blank in runs without explicit
insert boundaries, where the legacy whole-amplicon grade remains available.

## Demultiplex without consensus

Use the existing validation/run tools with `demux_only: true`; there is no
separate demux server to install. Alternatively save YAML with
`workflow: demux_only` and omit the flag on both tools. Reference libraries may
be omitted; barcode panels and input reads are still required for barcode sorting.
The CLI equivalent is `nanopore3 run --config my-run.yaml --demux-only`.

After updating the checkout, install into the same Python environment that your
AI host uses, then restart its MCP server:

```sh
python -m pip install -e ".[report,mcp]"
```

The server checks source identity; an already-running server must be restarted
when package code changes. Call `workspace_info` and read `nanopore3://guide`
after reconnecting. The input schemas for `start_validation` and `start_run`
should now include `demux_only`.

For a config at `my-run.yaml`, use these **tool calls** (JSON arguments, not shell
commands). Wait for each submitted job to reach `succeeded` before the next step:

```text
start_validation {"config_path":"my-run.yaml","quick":true,"demux_only":true}
job_status       {"job_id":"<validation job ID>"}
start_run        {"config_path":"my-run.yaml","run_id":"demux-preview","demux_only":true}
job_status       {"job_id":"<run job ID>"}
run_summary      {"run_dir":"runs/demux-preview"}
read_artifact    {"path":"runs/demux-preview/stages/02_demux/reads/index.csv"}
```

Poll `job_status` every few seconds; if it fails, report its error logs instead
of treating output presence as completion. `run_summary.workflow` is
`demux_only`, and a completed run has exactly `01_ingest` and `02_demux`.
There are no reference-assignment, consensus, protein-QC or culture-deconvolution
results. Read-length and mean-quality filters still apply.

| `run_summary.outputs` key | What the AI should inspect |
|---|---|
| `demux_fastq_index` | CSV of `scope`, `plate_id`, `well_id`, `reads`, `file`; FASTQ paths are relative to the index's directory |
| `demux_report` | HTML demux report |
| `demux_calls` | Gzipped call table, including failed/ambiguous calls and filter reasons |
| `provenance` | Configuration, input checksums, resource plan and source identity |

An index row with `file=by_well/P1/A1.fastq.gz` resolves to
`runs/demux-preview/stages/02_demux/reads/by_well/P1/A1.fastq.gz`. Preview it with:

```text
read_artifact {"path":"runs/demux-preview/stages/02_demux/reads/by_well/P1/A1.fastq.gz","max_chars":2000}
```

Use the paths actually returned by the index: identifiers may need filename
sanitisation and unoccupied bins have no FASTQ. `read_artifact` transparently
reads gzip text but returns a bounded preview; use the index's counts rather
than counting records in a truncated preview. All paths refer to the **server's
workspace**, not the AI host's filesystem.

Plate files include unresolved wells and **overlap** with well/unresolved files.
Never sum or concatenate all scopes as independent reads. To browse an entire
plate use its `plate` row; to browse one well use its `well` row. Report
`unresolved_well` counts separately. These are full oriented reads, not consensus
sequences; original IDs, sample IDs and unique read IDs remain in the headers.
Omit the well panel for plate-only output. See the
[configuration guide](configuration.md#demux-only-read-export) for filtering,
orientation and storage details.

To resume an interrupted demux run, use the same config/source/input and
`start_run(..., run_id="demux-preview", resume=True, demux_only=True)`. Do not
switch an existing demux run to full analysis with `resume`; choose a new run ID
and provide the reference configuration for a fresh full run from FASTQ.
