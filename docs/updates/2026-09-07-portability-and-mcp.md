# September 2026 update: portability, safe reuse and AI interfaces

This update adds an optional MCP server, repairs four reliability issues, and
aligns the CLI, local Jupyter and Colab workflows. Package version remains
`0.3.0` (pre-alpha); this is a repository update, not a PyPI release. Install the
updated checkout or a wheel built from it. The repository is now named
[PLACE-seq](https://github.com/t-j-fryer/PLACE-seq). When using private forks,
Colab needs an authenticated checkout or an uploaded wheel to use private code.

## Before updating an existing analysis

**Old results remain readable, but reuse requires the same source identity.**
Stage fingerprints now include the installed Python source and bundled presets,
including local edits. Runs without this identity, or made with different code,
cannot be resumed/inherited by this implementation. Recompute from the original
inputs into a new run ID. Do not change old manifests to bypass the checks.

Pin a commit or retain the exact wheel and dependency versions for future reuse.
Restart Python/notebook kernels and MCP servers after changing installed code.
Keep complete run directories for reuse; a report/summary export is insufficient.

## What changed

| Area | Behavior now |
|---|---|
| Subsampling | Validates the requested FASTQ prefix, rejects existing outputs and input aliases, and publishes atomically without clobbering another writer |
| Source provenance | Records installed source identity in run metadata/manifests and checks it before resume/rerun |
| Zero read cap | `consensus.maximum_reads: 0` includes every eligible read, including chimera members; positive caps retain their prior selection rules |
| CPU allocation | Defaults to process workers on all allocated CPUs; respects exposed affinity and Linux cgroup quotas/cpusets and clamps native threads |
| Group memory | Groups reads and output accumulators on disk; oversized retained groups fail clearly instead of silently dropping reads |
| CLI | Help/diagnostics expose source compatibility, effective CPUs and configured/default memory budgets |
| Notebook/Colab | Uses the kernel's interpreter, stops cells on CLI failures, creates fresh paths, preserves gzip inputs and archives complete reports with linked figures |
| MCP | 13 typed tools for project/configuration access, safe subsampling, validation, runs, reruns, polling, cancellation and result inspection |

Maximum parallelism remains the default:

```yaml
parallel:
  backend: process
  jobs: 0
  threads_per_job: 1
  group_memory_mb: 512
```

`group_memory_mb` guards **estimated retained data per group**, not total RSS.
Visible Linux cgroup memory limits may lower that budget. Reference indexes,
worker batches, native aligners and reporting still need additional memory.
See the [configuration guide](../configuration.md) for the estimate and controls.

Subsample publication needs hard-link support in the output directory. Use local
disk (such as APFS/ext4/NTFS); in Colab, stage on the VM disk before copying to
Drive. Filesystems that lack this operation fail safely.

## Install and run

From the updated checkout, using Python 3.10–3.13:

```sh
python -m pip install -e ".[report,mcp]"
python -m nanopore3 doctor
python -m nanopore3 init nanopore3-example
python -m nanopore3 run --config nanopore3-example/configs/example.yaml
```

For local Jupyter, add the `notebook` extra. The
[Jupyter/Colab notebook](../../notebooks/Nanopore3_Colab.ipynb) prefers a nearby
checkout during Setup and accepts a wheel/pinned package specification. Its MCP
cell uses that same installation rather than reinstalling a branch mid-analysis.

To let an MCP-capable AI host launch the server:

```sh
python -m nanopore3.mcp_server --workspace /absolute/path/to/project
```

The workspace must exist. Use the interpreter from the environment where MCP was
installed; Windows users can pass a Windows path as one quoted argument. For
client configuration, HTTP, Windows, Colab, tools and troubleshooting, see the
[MCP guide](../mcp.md). A chat window without MCP support needs a compatible host.
This update does not deploy a public server or configure an AI account.

## Evidence and limitations

Before publication, all **443 local tests passed** on macOS/Python 3.12.10,
including notebook execution and real stdio/HTTP MCP sessions. The installed
wheel also completed the MCP smoke workflow outside the checkout. Lint,
compilation and diff checks passed. GitHub Actions runs the portable suite on
Linux (Python 3.10/3.13), macOS and Windows, plus a Linux installed-wheel workflow;
[pull request #1](https://github.com/t-j-fryer/PLACE-seq/pull/1) records the
publication and checks. Live Colab is not covered by that matrix; its Drive mount
branch is simulated in tests.

The [hosted run for implementation commit 827da1e](https://github.com/t-j-fryer/PLACE-seq/actions/runs/34153711406)
completed successfully on every target:

| Hosted target | Result |
|---|---|
| Ubuntu / Python 3.10 | Passed: portable suite, MCP tests and wheel build |
| Ubuntu / Python 3.13 | Passed: lint, portable suite, MCP tests and wheel build |
| macOS / Python 3.12 | Passed: portable suite, MCP tests and wheel build |
| Windows / Python 3.12 | Passed: portable suite, MCP tests and wheel build |
| Linux installed wheel / Python 3.12 | Passed: reporting/MCP extras and real stdio workflow outside the checkout |

This table records that specific source revision; subsequent pull-request/main
checks remain available in GitHub Actions.

The final core-fix benchmark ran 120,000 synthetic reads with 16 automatically
selected workers. All nine checked artifacts matched the earlier serial results,
including contributor IDs, consensus summary counts, FASTA and QC. Median elapsed
time was 7.429 seconds, versus the historical 6.981-second auto-worker median.
Those measurements are unpaired and are not a production throughput forecast.

Read the [repository assessment](../repository-assessment.md),
[raw verification data](../benchmarks/2026-09-07-priority-fixes.json) and
[lab notebook](../../LAB_NOTEBOOK.md) for measurements, commands, compatibility
decisions and remaining work. No experimental reads or generated run outputs are
included in this update. Separate local RP04 clonality edits are excluded.
