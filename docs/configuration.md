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

### Clonality: how many clones a well should hold

A well may legitimately hold many clones. Where the design put **one clone per
block**, the expected diversity follows from the layout and is derived
automatically; where colonies were scraped, no expectation exists and inventing
one would turn an unknown into a spurious deviation.

```yaml
compressed_pcr:
  clonality:
    RP06: unspecified   # scraped with a multichannel: variable, report observed only
    RP07: per_block     # monoclonal picks mixed per block: expected = blocks pooled
```

`per_block` sums the blocks across that plate's pooled culture plates. Three
culture plates holding blocks 1-3, 4-6 and 7-10 give an expected 3 + 3 + 4 = 10
distinct genes per well, drawn on the clonality figure as a reference line.

Block keys are **nested by reference library**, because block numbering restarts
in each one — block 1 of `sumo_lab` is unrelated to block 1 of `aaseq_biotin`:

```yaml
  blocks:
    sumo_lab:
      "1": [CP_SUMO_A]
      "4": [CP_SUMO_B]
```

## Figures

Every run writes publication figures to `stages/06_report/figures/` as vector PDF
plus 600 dpi PNG, sized for a print journal (89 mm single / 183 mm double column,
6-8 pt sans). They are best-effort: a run does not fail when matplotlib is
absent, and `figures.json` records what was written or why it was skipped.

| Figure | Shows |
| --- | --- |
| `fig1_plate_occupancy` | distinct genes recovered per well, one 96-well panel per plate |
| `fig2_clonality` | observed clones per well by plate, against the expected count |
| `fig3_deconvolution` | read fate per plate; only when compressed PCR is enabled |

Regenerate without re-running the pipeline:

```bash
python scripts/make_figures.py --run runs/<run-id> --config configs/runs/<profile>.yaml
```

The categorical palette is **validated, not chosen by eye**. `#0072B2`,
`#D55E00`, `#009E73` clear the all-pairs colour-vision-deficiency check (worst
OKLab dE 11.0 under simulated protanopia and deuteranopia, 18.7 unsimulated). A
fourth hue drops the worst pair to 7.6, below target, so three is the cap wherever
two marks can touch. Re-check any change:

```bash
python scripts/validate_palette.py "#0072B2,#D55E00,#009E73" --pairs all
```

## Graded consensus output

Every run writes a browsable copy of the consensuses to
`stages/05_qc/consensus_by_plate/`, one FASTA per consensus:

```text
RP06/A01/RP06_A01__Block_1_dTF083_156_2__perfect.fasta
RP06/A01/RP06_A01__Block_7_dTF083_152_4__screenable.fasta
RP06/A01/RP06_A01__Block_9_dTF083_282_1__frameshift.fasta
```

Filed under `<plate barcode>/<well>/`, named `<plate>_<well>__<design>__<grade>`,
so a screening decision can be made from the listing without opening a file. The
name is self-describing if a file is moved. Wells are zero-padded so `A01` sorts
before `A10`, and a well holding several designs simply gets several files —
polyclonal wells are the normal case, not an error.

It lives under `05_qc` rather than `04_consensus` because the grade needs QC
results, and a promoted stage directory is immutable.

The grade collapses the QC criteria into one word. Precedence runs from "no usable
data" through "the construct is broken" to "the construct is fine", so the grade
reported is the most actionable problem rather than the first one found:

A well holding several designs is the **normal** case and simply produces several
files; `mixed_variants` is per design, not per well.

| Grade | Meaning |
| --- | --- |
| `perfect` | exact match to the designed reference, in frame, no internal stop |
| `screenable` | full length and in frame with no internal stop, but carries substitutions |
| `mixed_variants` | reads for this one design disagree beyond the support threshold: not a single clean clone |
| `mismatched` | identity or coverage below the QC floor |
| `truncated` | length outside the configured tolerance |
| `premature_stop` | a stop codon before the end of the reading frame |
| `frameshift` | the assembled reading frame is not a whole number of codons |
| `low_depth` | fewer contributing reads than `consensus.minimum_depth` |

Alongside the tree, `index.csv` lists every consensus with its grade, design,
culture plate, depth, identity and file path — including those with **no** file,
so a missing FASTA never has to be read as an oversight. `summary.json` gives
grade counts overall and per plate.

Rebuild without re-running the pipeline:

```bash
python scripts/export_consensus_tree.py --run runs/<run-id>
```

### Culture-plate recovery

Where compressed-PCR deconvolution is configured, `fig4_culture_plates` shows what
the compression actually recovered, and the report gains a **Culture plate
recovery** section with the pooled count, plates recovered, clones, median source
plates per well, and the weakest plate.

Each pooled colony-PCR plate gets two panels:

- a 96-well map of **how many distinct source culture plates each well recovered**,
  which answers "did this well see everything pooled into it";
- a bar per source culture plate of **clones recovered**, which answers "did every
  pooled plate contribute". A plate contributing far below its peers is coloured
  and labelled directly, since colour alone is not an encoding.

The pooled denominator comes from `compressed_pcr.pcr_plates`, not from what was
observed: a culture plate that contributed nothing must show as missing rather
than silently shrinking the denominator.
