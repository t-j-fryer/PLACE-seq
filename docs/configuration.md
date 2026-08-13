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

## Compressed PCR: gene identity as an extra demultiplexing key

When colonies from several culture plates are pooled into one colony-PCR plate,
the forward and reverse barcodes only identify a *PCR* well. That well may hold
colonies drawn from several culture plates at the same position. The assigned
gene supplies the missing coordinate, because a gene belongs to exactly one
assembly block and a block was picked into known culture plates:

```text
gene -> assembly block -> culture plate
```

Enable it by describing the pooling layout:

```yaml
compressed_pcr:
  enabled: true
  # Culture plates pooled into each colony PCR plate, keyed by plate barcode.
  pcr_plates:
    RP01: [CP_A, CP_B, CP_C]
    RP02: [CP_D]
  # Culture plate(s) each assembly block was picked into.
  blocks:
    "1": [CP_A]
    "2": [CP_B]
    "3": [CP_C, CP_D]      # split across two plates, in different PCR plates
```

A block split across culture plates is fine **provided those plates go to
different colony PCR plates**, so the reverse barcode separates them. If two
culture plates holding the same block are pooled into one PCR plate, no evidence
can tell them apart; that layout is rejected at configuration time rather than
producing ambiguous reads discovered after a run:

```text
error: compressed PCR layout is unresolvable: block 3 occupies CP_C, CP_D which
were both pooled into RP01. Split these culture plates across different colony
PCR plates so the reverse barcode separates them.
```

### Where the block comes from

By default the block is read from the reference identifier with
`block_pattern`, whose default `^Block_(\d+)_` matches identifiers such as
`Block_11_binder_dLK39_458_model`. Supplying the oPool design table per library
is authoritative and is joined **by sequence**, because the FASTA and the design
table use different identifiers:

```yaml
reference_libraries:
  aaseq_biotin:
    fasta: .../01_AAseq_Biotin_references.fasta
    fragments_csv: .../01_AAseq_Biotin_FULL_INFO.csv
```

### What it reports

`assignment_calls.csv.gz` gains `assembly_block`, `culture_plate`, and
`culture_plate_status`:

| Status | Meaning |
| --- | --- |
| `resolved` | the block names exactly one pooled culture plate |
| `unexpected_block` | the gene's block was never pooled into this PCR plate — a mis-pick, cross-contamination, or a misdescribed layout |
| `ambiguous` | several candidate plates; never attributed to one arbitrarily |
| `unknown_block` | the read was unassigned, or no block is known for its reference |
| `unknown_pcr_plate` | the plate barcode is not described by `pcr_plates` |
| `not_configured` | compressed-PCR mode is off |

`consensus.csv.gz` and the consensus FASTA carry `culture_plate` as provenance.
The consensus **grouping key is deliberately unchanged**, so consensus identities
are identical whether or not deconvolution runs; the source plate is recorded
alongside rather than folded into identity.
