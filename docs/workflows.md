# Worked examples

**Every run is amplicon data.** That never changes. Three things do:

| what varies | options | what it changes |
|---|---|---|
| **how much the primers span** | the expression cassette, or nearly the whole vector | only how much constant sequence flanks the designed region — the same mechanism, a different length |
| **what your reference list is** | the designed inserts alone, or sequences already assembled around them | how you describe where the designed region sits |
| **what a barcode holds** | one culture plate, or several pooled into one colony-PCR plate | whether the source plate has to be recovered from the gene |

These are independent. The four examples below are the corners:

| your references are… | one culture plate per barcode | several plates pooled |
|---|---|---|
| **the designed inserts** | [1](#1-inserts-one-plate-per-barcode) | [2](#2-inserts-pooled-colony-pcr) |
| **already assembled** | [3](#3-assembled-references-one-plate-per-barcode) | [4](#4-assembled-references-pooled-colony-pcr) |

Two things worth saying plainly, because the naming elsewhere has implied otherwise:

- **"Cassette" and "whole vector" are not modes.** They differ only in how many
  constant bases sit either side of the designed region. A 200 nt flank and an
  850 nt flank are described the same way and handled by the same code.
- **The designed region exists in both.** Whether your amplicon covers the cassette
  or the whole vector, the insert — your oligo-pool library, or any reference list
  you supply — is somewhere inside it, and is what gets scored separately.

Before any example:

```bash
python -m pip install -e .
nanopore3 doctor
```

---

## 1. Inserts, one plate per barcode

**You have:** a FASTA of designed inserts (an oligo pool order), one culture plate
per barcode, one colony per well.

**Decide one thing:** do you want the constant sequence around the insert
reconstructed and scored, or only the designed region?

### 1a. Score only the designed region

Say nothing about flanks. The consensus covers what the primers bound; the
reference is the insert, so that is what accuracy is measured over.

```yaml
schema_version: 1
preset: ont-r10-amplicon
run_name: cassette-run
output_root: ./runs

inputs:
  - path: ${SEQ_DATA}/run1.fastq
    sample_id: run1

reference_libraries:
  designs:
    fasta: references/oligo_pool_designs.fasta

plate_reference_map:
  BC01: designs

library:
  name: cassette-amplicon
  forward_motif: CATAATCCGCACGCATCTGG      # bounds the insert
  reverse_motif: CGGATTGGCGAATGGGACGC
  minimum_read_length: 300
  maximum_read_length: 1200

barcodes:
  plate:
    registry_csv: barcodes/plate_barcodes.csv
    family_id: my_plates
  well:
    registry_csv: barcodes/well_barcodes.csv
    family_id: my_wells

random_seed: 1
```

### 1b. Also reconstruct and score the constant regions

Add the constant context. This works identically whether that context is 200 nt of
cassette or 850 nt of vector:

```yaml
reference_libraries:
  designs:
    fasta: references/oligo_pool_designs.fasta
    flanks:
      template: references/example_construct.fasta   # one library member, already assembled
      anchor_length: 20

qc:
  # Locate the ORF inside the consensus: upstream_constant begins at the start
  # codon, and the tail of downstream_constant ends at the terminal stop.
  upstream_constant: ATGCAGCTT
  downstream_constant: ...GGATCCCATCACCACCACCATCACTAA
```

`flanks.upstream` / `flanks.downstream` are the alternative if you would rather
paste the two constant regions than supply an example construct.

**What you gain:** `insert_identity` and `insert_edit_distance` reported separately
from the flanks, acceptance floors measured over the insert rather than diluted by
constant sequence, and chimera detection scoped to the region that actually
discriminates.

```bash
export SEQ_DATA=/path/to/sequencing/data
nanopore3 validate --config run.yaml
nanopore3 run --config run.yaml
```

---

## 2. Inserts, pooled colony PCR

**You have:** the same reference list, but several culture plates were combined into
one 96-well colony-PCR plate, so a well coordinate no longer identifies a source
plate. What resolves it is the gene: a design only ever picked onto plate 3 came
from plate 3.

Everything from example 1 (either variant — pooling is independent of whether you
reconstruct the flanks), plus:

```yaml
compressed_pcr:
  enabled: true
  block_pattern: "^(Block_\\d+)_"     # extracts the block from a reference name
  layout_csv: my_layout.csv
```

`my_layout.csv` — one row per (culture plate, block), read as *this block went on
this plate, and this plate went into this barcode*:

```csv
plate_barcode,culture_plate,library,block
BC01,PlateA,designs,Block_1
BC01,PlateB,designs,Block_2
BC01,PlateB,designs,Block_3
```

This table is **experimental design**. You supply it; it is never inferred from the
reads, because inferring it would use the data to build the key that then
interprets that same data.

Check it before spending an hour:

```bash
nanopore3 layout --config run.yaml
```

which prints the layout and names any block with no plate, plate with no block, or
barcode missing from `plate_reference_map`.

**What changes:** clones are filed under `<barcode>/<culture plate>/<well>/`, and
`index.csv` gains a populated `culture_plate`. A well that cannot be attributed to
one plate keeps the `|`-joined candidates and stays at well level — unresolved is
reported as unresolved, never guessed.

---

## 3. Assembled references, one plate per barcode

**You have:** the library already assembled — every expression cassette, or every
whole vector, as a complete sequence. **These are the same case.** The pipeline
finds the shared backbone; whether it is 200 nt or 3 kb makes no difference.

```yaml
reference_libraries:
  designs:
    fasta: references/assembled_library.fasta
    flanks:
      derive: true
      anchor_length: 20
```

`derive` takes the longest common prefix and suffix across the set as the constant
regions; what remains is the designed region. No alignment is involved — members of
one library share exact sequence at both ends.

It refuses rather than guesses when there are fewer than two references, under
60 nt shared at either end, or no variable region between the shared ends.

**Or skip `flanks` entirely** and score each construct as a whole. You lose
per-region accuracy, insert-scoped floors and insert-scoped chimera detection — see
[References](references.md#why-this-matters-beyond-tidiness) for why that is worth
avoiding.

Everything else is exactly as example 1.

---

## 4. Assembled references, pooled colony PCR

Example 3 plus the `compressed_pcr` block from example 2. Nothing interacts: the
layout table resolves the source plate, `derive` locates the designed region, and
neither knows about the other.

The shipped version of this is
[`configs/runs/260608_full_length_short.yaml`](../configs/runs/260608_full_length_short.yaml)
— though it supplies inserts plus a template rather than assembled references, so it
is really example 2 with a whole-vector amplicon.

---

## Reading the results

`consensus_by_plate/index.csv` is where most analysis should start. One row per
clone:

| column | use |
|---|---|
| `grade` | Legacy whole-amplicon grade: `perfect`, `screenable`, `mixed_damage`, … |
| `insert_grade`, `vector_status` | Separate bounded-insert and sequenced-vector results when [explicit QC boundaries](configuration.md#separate-insert-and-vector-qc) are configured |
| `culture_plate`, `well_id` | where to go back to on the bench |
| `insert_identity`, `insert_edit_distance` | accuracy of the designed region alone |
| `flank_5p_edit_distance`, `flank_3p_edit_distance` | errors in the constant sequence |
| `designed_allele_fraction` | for a mixed clone: how much of the well still carries the design |
| `mixed_worst_effect` | `silent` / `missense` / `nonsense`, or empty outside the ORF |
| `file` | the FASTA, relative to the tree root |

A `mixed_damage` clone is a real clone carrying a pre-transformation lesion, not an
artefact: usable, with the caveat those last two columns spell out.

---

## Changing your mind afterwards

Retuning a threshold does not mean repeating the demultiplexing, and does not need
the FASTQ still attached:

```bash
nanopore3 rerun --config tuned.yaml --from-run runs/<run-id> --from 04_consensus
```

Stages before `--from` are carried into a new run; everything after is recomputed.
A change that reaches back into an inherited stage is refused by name rather than
silently built upon. See [Architecture](architecture.md#reruns-recomputing-the-analysis-without-the-fastq).

## Inspect demultiplexed reads without making consensus

Use `python -m nanopore3 run --config my-run.yaml --demux-only --run-id demux-preview`.
Open `stages/02_demux/report.html` and `stages/02_demux/reads/index.csv` in that run.
The index links plate and well FASTQ bins; references are optional. For a notebook
or Colab, set `DEMUX_ONLY = True`. For MCP, pass `demux_only=True` to validation and
run tools. See [demux-only export](configuration.md#demux-only-read-export) for
filtering, orientation, unresolved wells and overlapping output views.
