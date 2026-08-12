# Configuration guide

Nanopore3 resolves every path relative to the YAML configuration file and rejects unknown keys.
This makes configurations portable while preventing a misspelled threshold from silently using a
default. Start with `nanopore3 init` and edit the generated example.

## Main sections

- `inputs`: one or more FASTQ/FASTQ.gz paths with stable sample IDs.
- `references`: one legacy/default FASTA set, or `reference_libraries` plus a
  `plate_reference_map` for multiplexed experiments.
- `library`: required ordered boundary motifs and optional length/mean-quality filters.
- `barcodes.plate` and `barcodes.well`: inline sequences or one CSV
  `registry_csv`/`family_id`, plus search ends, edit limits and an explicit
  decision policy.
- `parallel`: `auto`, `thread`, `process`, or `serial`, with a bounded job,
  chunk-size and nested-thread budget.
- `consensus`: explicit `portable` or optional `mafft_spoa` backend, depth cap,
  minimum depth and support fraction.
- `qc`: independent whole-sequence identity, coverage and length thresholds.

Configured motifs are mandatory evidence. A missing pair produces `motif_missing`; it does not
fall back to full-read/vector assignment. References with identical sequences are emitted as one
alias set because the sequence evidence cannot distinguish their names.

With more than one reference library, every processed plate must be mapped. An
unmapped plate remains visible as `unmapped_reference_library` and is excluded
from consensus; reference sets are never pooled as an implicit fallback. CSV
registry paths are resolved relative to the YAML, and their contents and digest
are included in the resolved run configuration.

## Conservative assignment

The reference index counts distinct canonical query k-mers and weights each k-mer by the number
of reference-sequence groups that own it. It shortlists candidates, then validates them with
semi-global edlib alignment. A call must pass identity, query coverage, reference coverage and a
best-versus-second identity margin. Calls that fail remain explicit in the per-read table and
never enter consensus construction.

Smaller `kmer_sizes` are fallback modes and are tried in the listed order. Thresholds should be
selected on labeled or synthetic development data, then evaluated once on held-out runs. Do not
tune thresholds using the same experimental run being reported.

## Parallelism

`parallel.jobs` controls concurrent read batches and `chunk_reads` controls the
number of records in each batch. Barcode panels are compiled once, batches are
bounded and ordered, the parent owns output writing, and result order is
independent of task completion. `thread` avoids notebook/process startup issues;
explicit `process` uses spawn-safe module workers for CPU scaling from the CLI.
Use `serial` when diagnosing a restricted notebook environment. Do not assume
that every logical core is fastest: benchmark the fixed assay profile on the
target machine and select the throughput knee.

The test suite verifies that serial and threaded runs produce identical normalized scientific
artifacts. Throughput and memory still need benchmarking on representative full-sized data before
the pre-alpha label is removed.
