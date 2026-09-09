# Configuration guide

PLACE-seq resolves every path relative to the YAML configuration file and rejects unknown keys.
This makes configurations portable while preventing a misspelled threshold from silently using a
default. Start with `nanopore3 init` and edit the generated example.

## Running on a machine you did not configure for

Parallel execution defaults to `backend: process`, `jobs: 0` and
`threads_per_job: 1`: use all detected CPUs for the parallel stages. An explicit
positive `jobs` value still caps the worker count, and `backend: serial` disables
process parallelism. On a single-CPU machine, execution falls back to the local
worker without requiring process initialization.

`parallel.jobs: 0` uses the effective CPU allocation. Detection considers the
host count, process CPU count/affinity when the OS exposes it, and Linux cgroup
v1/v2 CPU quotas and cpusets, including visible parent limits. Fractional quotas
round down to a whole worker, with a minimum of one. Both workers and native
threads are clamped so `jobs * threads_per_job` stays within this budget. The
bounded thread count is passed to MAFFT. `doctor --json` reports both host
`cpu_count` and effective `available_cpus`, plus `parallel_defaults`,
`group_memory_limit_bytes` and `analysis_implementation`. These are default
settings; `validate --config ... --quick` reports the configured CPU plan under
`resources` and the configured/resolved group budget under
`group_memory_limit_bytes`. MCP exposes the same doctor fields through
`workspace_info.runtime` and the same validation JSON through job status. Linux hierarchy details are described
in the [kernel cgroup v2 documentation](https://docs.kernel.org/admin-guide/cgroup-v2.html).

Consensus read groups, chimera wells, claimed IDs and their output accumulators
use temporary SQLite storage under the active stage. Only one group is processed
at a time; read order and deterministic selection are preserved. Allow free local
disk space for uncompressed intermediates and the temporary database, especially
in Colab. Temporary databases are removed on success and handled failures.

`parallel.group_memory_mb` defaults to **512 MiB**. It limits the *estimated
retained working data for one group*, using 2,048 bytes plus 64 times the sequence
length per retained read. A visible Linux cgroup memory limit additionally caps
this budget at one quarter of that limit. An oversized group fails explicitly
before stage publication; reads are never silently discarded to fit memory.
This is admission control, **not a hard limit on total RSS**: reference indexes,
worker copies/batches, native aligners and later reporting also consume RAM. For
those costs use a smaller explicit worker/chunk count when needed. OS memory
limits remain the mechanism for enforcing total process/container memory.

`consensus.maximum_reads: 0` uses every eligible read, including chimera clone
members. A positive value retains the existing deterministic cap. Unlimited
groups can hit the memory guard; raise `group_memory_mb` only with sufficient RAM
or split the analysis. Resource settings change execution, not scientific filters.

Progress is reported per stage by default — a run that prints nothing for an hour is
indistinguishable from one that has hung, and on a hosted notebook it also risks
being disconnected for idleness. `--quiet` suppresses it.

`nanopore3 subsample --input X --output Y --reads 100000` takes the first N reads of a
FASTQ, for trying a configuration or checking barcode recovery before committing to a
full run. Reads come in file order, so it is reproducible and costs one pass over the
head of the file rather than over all of it. The output must be a new path:
existing files and input aliases are refused. A validated prefix is published
atomically without overwriting a competing writer's output. This requires hard
link support in the output directory (e.g. local APFS/ext4/NTFS); unsupported
filesystems fail safely. In Colab, write the subset on the local VM disk before
copying it to Drive.

## Paths that travel

A configuration describes an experiment. *Where* that experiment's data sits is a
property of the machine reading it, and writing the second into the first is what
makes a config unshareable — and, once committed, what puts a home directory in a
public repository.

Use `${NAME}` for anything machine-specific:

```yaml
inputs:
  - path: ${SEQ_DATA}/260608/AI_DBTL.fastq
```

```bash
export SEQ_DATA=/Volumes/MyDrive
```

Values come from the environment. An unset variable is an error naming it, never an
empty string — expanding to `""` would give `/AI_DBTL.fastq` and a file-not-found
that says nothing about the cause:

```
error: inputs[0].path refers to unset environment variable(s): SEQ_DATA.
Set them for this machine, for example:
  export SEQ_DATA=/path/to/data
```

Relative paths resolve against the **configuration file**, not the working
directory, so a config next to its references needs no variables at all.

Anything genuinely private — unpublished designs, collaborator data — belongs in
`configs/local/`, which is git-ignored.

## Presets: the tuning you cannot set from first principles

Roughly a quarter of a real configuration is platform tuning: k-mer sizes,
significance, minor-allele floors, chimera windows. A new user has no basis for any
of it. `preset:` supplies a named, versioned bundle:

```yaml
preset: ont-r10-amplicon
```

There is one: `ont-r10-amplicon`. The input is always an amplicon, and neither the
span of the primers nor the form of the reference list changes the tuning — so a
second preset would only have implied a distinction that does not exist.

The preset supplies tuning; your configuration supplies the experiment and **wins
wherever both mention a key** — including one key inside a section, leaving the rest
of that section inherited:

```yaml
preset: ont-r10-amplicon
consensus:
  maximum_reads: 50      # this wins; minimum_depth, significance etc. still inherited
```

A preset is a way of *writing* a configuration, never a hidden layer under one. The
merged result is what is validated, digested and recorded in the run, so two runs
naming the same preset are exactly as reproducible as two spelling it out.

`configs/runs/260608_full_length_short.yaml` is the shipped 260608 run written this
way: **222 lines to 117**, resolving to identical settings.

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

## Graded clone output

`consensus_tree: true` (the default) writes `<run>/consensus_by_plate/`: one FASTA
per clone, filed under `<barcode>/<culture plate>/<well>/` and named with its design
and a single-word grade, plus `index.csv` — every clone, one row each, with grade,
identity, per-region accuracy and the mixed-allele annotation. That table is what
most downstream screening analysis should read.

Set it to `false` only for a very large run where writing one small file per clone
is itself a cost; the same tree can then be built afterwards with
`python scripts/export_consensus_tree.py --run <run>`.

The tree is **regenerated on every run and rerun**, so do not keep your own files
inside it — a stale FASTA is indistinguishable from a current one, which is why it is
replaced rather than merged. A directory the exporter did not write (one with no
`index.csv`) is refused rather than replaced, so your own directory of that name is
safe.

## Compressed PCR: gene identity as an extra demultiplexing key

> The layout can be supplied as a spreadsheet rather than nested YAML, and is
> cross-checked against the references at preflight. See
> **[Pooling layout](pooling-layout.md)** for the table format, the monoclonal
> case (which needs no configuration at all), and the migration command.

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
consensus sequences per well, drawn on the clonality figure as a reference line.

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
| `fig1_plate_occupancy` | **consensus sequences** per well, one 96-well panel per plate |
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

## Separate insert and vector QC

For full-vector amplicons, define the insert explicitly with a pair of unique,
forward-oriented reference motifs. The insert is **strictly between** these
motifs; both boundary motifs belong to the vector. Each motif must occur once
in every analysed reference, in the given order. Validation rejects missing,
repeated or reversed boundaries.

```yaml
qc:
  insert_left_boundary: CCGATGCAGCTT  # synthetic example; replace for your assay
  insert_right_boundary: TCGGGCCACTAA
  upstream_constant: ATGCAGCTT       # begins at the actual ORF start codon
  downstream_constant: TCGGGCCACTAA  # includes the terminal stop
```

Use this same YAML through the CLI, local notebook, Colab or MCP `save_config` /
`start_validation` / `start_run` tools. No separate interface option is needed.
The motifs must remain inside the analysed reference after extraction-anchor
trimming. Without explicit insert boundaries, existing grading and filenames are
retained; derived shared flanks alone do not enable these separate grades.

For reference-guided consensuses, the clone index and detailed QC table add:

| Field | Meaning |
| --- | --- |
| `insert_grade` | `perfect` means exact, unambiguous insert DNA. `screenable` means differing DNA that passes insert identity, coverage, length and coding checks. Other calls are `mixed_variants`, `frameshift`, `premature_stop`, `truncated`, `mismatched`, `low_depth` or `not_evaluable`. |
| `vector_status` | `perfect` means exact, unambiguous sequence outside the insert; `edited` means one or more differences or ambiguous bases; `not_evaluable` means unavailable sequence. |
| `insert_coding_status` | `pass`, `fail` or `not_evaluable`. An exact DNA insert does not establish protein function; without a configured reference reading frame its coding status is not evaluable. |
| `insert_edit_distance`, `vector_edit_distance` | Edit counts in each scope; vector includes both sequenced flanks. |
| `insert_ambiguous_bases`, `vector_ambiguous_bases` | Ambiguous consensus bases in each scope. |

FASTA names become `...__insert-perfect__vector-edited.fasta`. Both labels and
insert coding status also appear in headers. `summary.json` adds `insert_grades`
and `vector_statuses`; the HTML report presents both. Chimeric clones retain
their existing scaffold-based grade and do not receive these separate labels. The existing `grade` column
and `grades` summary remain whole-amplicon classifications for compatibility.
Whole-ORF `reading_frame` and `internal_stops` remain separate: an edit to a
boundary or tag can affect the ORF even if the bounded insert is perfect.

When both query boundary motifs are intact and unique, alignments are pinned to
them to prevent repetitive sequence from moving an edit across a boundary.
Otherwise the global alignment projects the reference boundaries onto the query.
An insertion immediately before the right intact motif belongs to the insert.
Insert coding checks use the reference ORF phase, supplementing partial edge
codons with reference bases; they describe the insert independently of backbone
changes. Ambiguous inserts are conservatively `mixed_variants`, without assigning
a biochemical cause. Neither a grade nor consensus depth proves sample purity.

Vector status covers the **sequenced amplicon outside the insert**, not unsequenced
plasmid regions or trimmed extraction anchors. For plate-specific backbone
variants, supply separate reference libraries and route the appropriate plate
with `plate_reference_map`. Keep original references and previous runs, record
the intended change, and run again from raw reads after changes to pipeline code.

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

## Source identity and safe reuse

Stage fingerprints and `run.json` now include a SHA-256 identity of all installed
package Python source files and bundled presets, with relative paths and
normalized newlines. It works in a checkout or installed wheel without Git and
includes local edits. Documentation and bytecode cache changes do not affect it.
A long-running Python process must restart after package source changes.

Resume and rerun require the same implementation identity. Older runs lacking
it, or runs made with different code, remain readable but cannot be inherited or
resumed by this version. Start a new run with a new ID from the original inputs;
do not edit old manifests to bypass the check. Configuration-tuning reruns still
work with unchanged code when the inherited stage fingerprints match. This
conservative policy invalidates reuse even for code changes unrelated to a
particular stage; dependency/native-backend versions remain separately recorded.

## Demux-only read export

Use this when you want to inspect or analyse reads by barcode plate and well
without reference assignment, consensus building, protein QC or pooling
deconvolution. It runs only `01_ingest` and `02_demux`.

```bash
python -m nanopore3 validate --config my-run.yaml --demux-only --quick
python -m nanopore3 run --config my-run.yaml --demux-only --run-id demux-preview
```

Alternatively set `workflow: demux_only` in YAML and omit both CLI flags. The
default is `workflow: full`. The flag overrides YAML for that invocation and is
recorded in `run.json`. It must also be supplied on `--resume` if your YAML still
says `full`. Resume requires matching code, config and input checksums, and
verifies the FASTQ exports. Switching workflow is a config change: choose a
fresh run ID and run from FASTQ to perform a subsequent full analysis. The
`rerun` command is for full analyses, not demux-only workflows.

Minimal synthetic configuration (replace barcode sequences and filters with
your experiment's values):

```yaml
schema_version: 1
workflow: demux_only
output_root: runs
inputs:
  - path: reads.fastq.gz
    sample_id: sample1
library:
  minimum_read_length: 100
  maximum_read_length: 10000
  minimum_mean_quality: 10
barcodes:
  plate:
    sequences:
      P1: AAAACCCC
    max_edits: 0
    search_window: 30
    search_ends: [head]
    allow_reverse_complement: true
  well:
    sequences:
      A1: ACGTACGT
    max_edits: 0
    search_window: 60
    search_ends: [head]
    allow_reverse_complement: false
parallel:
  jobs: 0
```

References and plate-to-reference mappings may be omitted. An existing full
configuration also works; reference FASTAs and optional consensus executables
are not required or opened in this mode. YAML schema and barcode validation
still apply. Read-length and mean-quality gates remain active, including any
preset defaults. Omit `barcodes.well` for plate-only output; no placeholder well
FASTQs are then written. Without a plate panel, sample IDs provide the grouping.

```text
runs/demux-preview/
  run.json
  stages/
    01_ingest/
    02_demux/
      report.html
      summary.json
      demux_calls.csv.gz
      demuxed_reads.jsonl.gz
      reads/
        index.csv
        by_plate/P1.fastq.gz
        by_well/P1/A1.fastq.gz
        unresolved_well/P1.fastq.gz
```

| Output | Contents |
| --- | --- |
| `by_plate` | Reads passing length/quality gates with an accepted plate call, including unresolved or conflicting well calls. |
| `by_well` | Reads with accepted plate and well calls and consistent orientation. |
| `unresolved_well` | The subset of plate reads whose well call was not accepted. Created only when needed. |
| `index.csv` | Scope, original plate/well IDs, read count and FASTQ path relative to `reads/`. |
| `demux_calls.csv.gz` | One row per input read, including rejected/ambiguous calls and filter reasons. |
| `demuxed_reads.jsonl.gz` | Internal fully accepted-read stream; unresolved wells are excluded. |

The plate and well/unresolved files are **overlapping views of the same reads**;
do not combine them as independent inputs. Reads from multiple input files with
the same barcode IDs are pooled in the same output file. FASTQ headers preserve
the original read ID and add unique `read_uid` and `sample_id` provenance.
Sequences and quality scores are retained in full, with reverse complementation
and quality reversal when barcode orientation requires it. Primers/adapters are
not trimmed and reference identity is not checked. Unassigned plate reads and
reads failing length/quality filters remain accounted for in the call table but
are not exported to plate/well FASTQ files.

Only occupied bins produce FASTQs; a run with zero accepted reads still writes
an index header, call table and report. File components that need sanitisation,
case-collision protection or Windows reserved-name protection get a stable hash
suffix; consult the index for exact paths. Exporting uses at most 32 open FASTQ
handles and streams records rather than retaining all reads in memory. Files
reopened after eviction use concatenated gzip members supported by standard gzip
readers. Output files are checksummed stage artifacts and only published when the
whole demux stage succeeds.

In the **local notebook or Colab**, set `DEMUX_ONLY = True` in the configuration
cell. References are removed from the generated config. The result cell displays
the FASTQ index and the archive cell copies the demux outputs to durable storage.
For **MCP**, pass `demux_only: true` to `start_validation` and `start_run`, or use
YAML `workflow: demux_only`. Poll the job, then read `run_summary.outputs.demux_report`,
`demux_fastq_index` and `demux_calls`. A completed demux-only run has exactly two
stages; missing consensus and QC stages are intentional.
