# Nanopore3 architecture

## Purpose

Nanopore3 separates the useful algorithms developed in Nanopore2 from notebook state, hardcoded
paths, accumulated output directories, and assay-specific plots. The production path is a small
Python library with explicit stage contracts. A CLI, notebook, or future workflow engine invokes
the same functions and receives the same structured results.

This document describes the target contract. Features not yet exposed by `nanopore3 --help` are
planned v0.1 work, not claims about the current implementation.

## Pipeline model

```text
validated configuration and input manifest
                  |
                  v
plate demultiplex -> well demultiplex -> candidate assignment
                                           |
                                           v
                               alignment validation / rejection
                                           |
                                           v
                                      consensus
                                           |
                                           v
                               sequence QC -> reports
```

Each stage consumes declared inputs and emits a machine-readable table plus any sequence files.
The next stage reads the table rather than reconstructing biological meaning from directory or
file names.

The assignment vocabulary is explicit:

- `assigned`: sufficient identity, coverage, support, and score margin;
- `ambiguous`: plausible candidates cannot be separated;
- `motif_failed`: required boundaries cannot be located;
- `truncated`: evidence supports an incomplete expected construct;
- `chimera`: spatially separated, reference-specific evidence supports a breakpoint;
- `no_match`: no candidate passes the configured evidence floor.

Quality control likewise distinguishes `pass`, `fail`, and `not_evaluable`.

## Package boundaries

The intended modules have narrow responsibilities:

- configuration: schema validation, profiles, paths, and normalized hashes;
- I/O: strict FASTA/FASTQ streaming, sequence identifiers, and atomic files;
- demultiplexing: barcode evidence and orientation without report generation;
- assignment: candidate generation, validation, margins, and reason codes;
- consensus: deterministic read selection and pluggable consensus builders;
- QC: region, coding, reference, truncation, and chimera evidence;
- reporting: tables and visualizations derived only from canonical results;
- runtime: resource planning, backend discovery, provenance, and stage state.

Algorithms should be usable independently through Python. CLI functions translate arguments and
format errors; they do not contain the algorithms.

## Portability and backend policy

The mandatory backend uses pip-installable Python plus edlib and works on macOS, native Windows,
Linux, and Colab. Optional adapters may use:

- minimap2 or mappy for high-throughput candidate alignment;
- SPOA for partial-order consensus;
- MAFFT followed by an internal deterministic column caller for MSA consensus.

Backend discovery happens before a stage starts. Executables are invoked with argument arrays,
never shell interpolation, in isolated temporary directories. Their path and version are captured.
A named unavailable backend is an error. Automatic selection is allowed for convenience but may
not be silent: the selected implementation is written to provenance.

Different consensus algorithms may produce different scientifically valid results. A canonical
study must therefore name and lock its backend; portability does not imply bit-identical output
across different algorithms.

## Parallel execution

Nanopore3 uses a central resource budget rather than independent thread choices in every stage.
The main controls are:

- `jobs`: concurrent samples, chunks, or consensus groups;
- `threads_per_job`: threads given to an external tool;
- `backend`: `auto`, `serial`, `thread`, or `process`.

The planner enforces `jobs * threads_per_job` within the available CPU budget and can additionally
cap workers using an estimated memory-per-job. Defaults should be conservative; users can opt in
to more parallelism.

Process workers are module-level callables with serializable inputs and spawn-safe behavior.
Nothing starts a pool during import. Work is chunked coarsely to limit Windows spawn overhead.
Workers do not append concurrently to a canonical output: they return records or write unique
temporary parts, and the parent performs a sorted, deterministic merge. Notebook execution may
prefer threads or serial execution where process startup is fragile.

External tools introduce nested parallelism, so the scheduler never gives every MAFFT process all
CPUs. Serial and parallel runs over the same input must be identical in golden tests.

## Immutable run layout

A run is a new directory, conventionally identified by a timestamp and normalized-configuration
hash:

```text
runs/<run-id>/
  manifest.json
  run.log
  events.jsonl
  01_plate_demux/
  02_well_demux/
  03_assignment/
  04_consensus/
  05_qc/
  06_report/
```

Inputs are read-only. Each stage writes a temporary directory and atomically promotes it only
after validation, then records completion in the manifest. Partial output is never interpreted as
a completed stage. Resume is permitted only when configuration and input checksums match. A new
configuration creates a new run rather than deleting or mixing previous results.

Use short, portable internal identifiers in paths. Full biological names belong in tables. This
avoids Windows path-length and forbidden-character problems and makes rename operations harmless.

## Provenance contract

The run manifest should contain:

- manifest schema and Nanopore3 version;
- source revision or source-tree hash;
- normalized configuration and its hash;
- Python, operating system, architecture, CPU budget, and package versions;
- selected backend names, executable paths, versions, and parameters;
- input/reference paths, sizes, and SHA-256 checksums;
- random seeds and deterministic sampling policy;
- stage start/end timestamps, status, record counts, and key output checksums.

Per-read assignment records preserve orientation, motif state, best and second candidates, their
scores, score margin, final state, and reason code. Consensus records preserve all contributor
read IDs, depth, backend, and support metrics. These records are the source of truth for reports.

CSV/CSV.gz and JSON are the portable interchange formats in the minimal installation. Parquet can
be offered with reporting dependencies but is not required to run the core pipeline.

## Configuration

Configuration separates assay facts from algorithm choices:

- an assay/library profile defines motifs, barcode sets, references, expected plate layout, and
  stable identifiers;
- a run file selects inputs, outputs, assay profile, thresholds, resources, and backends.

Validation occurs before processing. It checks sequence alphabets, duplicate identifiers,
barcode ambiguity, required motifs, file readability, output conflicts, and backend availability.
Defaults are visible in a resolved configuration stored with the run.

## Reproducibility levels

Nanopore3 supports three complementary environments:

1. Portable pip installation: easiest use and broadest platform support.
2. Platform lock: exact Python packages and native tools for macOS/Linux/WSL.
3. Versioned Linux container: canonical archival and high-throughput reproduction.

Package metadata alone does not lock native tools. A reported scientific run should cite its
manifest and, where used, environment lock or container digest.

## v0.1 scope and validation gate

Version 0.1 establishes configuration, strict sequence I/O, portable matching/alignment
primitives, deterministic runtime behavior, immutable run scaffolding, and synthetic tests. It
does not promise that every exploratory Nanopore2 branch has been migrated.

A classifier or consensus backend becomes recommended only after validation against synthetic
reads with known substitutions/indels/truncations/chimeras and held-out experimental controls.
The gate includes assignment precision/recall, false-chimera rate, consensus accuracy by depth,
runtime/memory measurements, and equivalent serial/parallel output.


## Reruns: recomputing the analysis without the FASTQ

Demultiplexing dominates a run's wall clock and depends only on the barcodes.
Changing a consensus threshold and repeating the whole thing wastes it, and
requires the original FASTQ to still be attached — which, for data on removable
media, it often is not.

```
nanopore3 rerun --config tuned.yaml --from-run runs/my-run --from 04_consensus
```

Stages before `--from` are carried into a **new** run directory and everything from
`--from` onward is recomputed. Two properties are preserved:

- **Runs stay immutable.** A rerun never edits its source. Inherited artifacts are
  hard linked where the filesystem allows, so carrying 600 MB of intermediates
  forward costs no disk and still cannot be modified in place.
- **Reuse is decided by the stage fingerprint, not by trust.** Every stage records
  a fingerprint over its parameters, input digests and pipeline version. A stage is
  inherited only if the current configuration reproduces that fingerprint. Change a
  barcode setting and try to rerun from `04_consensus`, and it is refused:

  ```
  error: the configuration differs from the one that produced an inherited stage.
  Inherited: 01_ingest, 02_demux, 03_assignment. Re-run with --from set to the
  earliest stage your change affects, or run from the FASTQ instead.
  ```

  The partially built run directory is removed, so a refusal leaves nothing behind.

The FASTQ is not required when both stages that read it (`01_ingest`, `02_demux`)
are inherited. Provenance still names the input file and its checksum: those are
carried from the source run's manifests rather than recomputed, and `run.json`
records `inherited` — the source run id, which stages came from it, and where
recomputation began.

A stage that did not run in the source (chimera detection, when it was disabled) is
not an error: it is computed if the new configuration calls for it.

### `--resume` versus `rerun`

`--resume` continues *the same run under the same configuration*, and refuses if
the configuration digest differs at all. Use it after an interruption. Use `rerun`
when the configuration has changed and you want the unaffected stages kept.
