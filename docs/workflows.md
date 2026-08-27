# Worked examples

Three complete scenarios, from what you have on the bench to what you read
afterwards. Each is a whole configuration — nothing is elided.

Before any of them:

```bash
python -m pip install -e .
nanopore3 doctor            # confirms the runtime and any optional backends
```

---

## 1. Monoclonal plates, insert-only amplicon

**You have:** one culture plate per barcode, one colony per well, and your primers
amplify the designed region itself. This is the simplest case and needs no
whole-vector handling and no deconvolution.

```yaml
schema_version: 1
preset: ont-r10-amplicon-insert
run_name: my-first-run
output_root: ./runs

inputs:
  - path: ${SEQ_DATA}/run1.fastq
    sample_id: run1

reference_libraries:
  designs:
    fasta: references/designs.fasta

plate_reference_map:
  BC01: designs

library:
  name: my-amplicon
  forward_motif: CATAATCCGCACGCATCTGG
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

```bash
export SEQ_DATA=/path/to/sequencing/data
nanopore3 validate --config my_run.yaml     # checks everything before spending time
nanopore3 run --config my_run.yaml
```

**You get** `runs/<run-id>/`:

```
consensus_by_plate/            one FASTA per clone, named with design and grade
  index.csv                    every clone, one row each  <- read this one
  BC01/A01/BC01_A01__design_17__perfect.fasta
stages/06_report/report.html   per-stage counts and QC
```

Filter `index.csv` on `grade == "perfect"` to get the wells worth keeping.

---

## 2. Monoclonal plates, whole vector

**You have:** the same, but your primers sit outside the insert, so each read
carries constant vector sequence at both ends. You want accuracy reported
separately for the designed part and the vector.

Everything above, plus a way to say where the insert is. Pick whichever matches
what you have (see [References](references.md) for all four):

```yaml
preset: ont-r10-amplicon        # note: not the -insert variant

reference_libraries:
  designs:
    fasta: references/designs.fasta      # inserts only
    flanks:
      template: references/example_construct.fasta   # one binder already in the vector
      anchor_length: 20

qc:
  # These locate the ORF inside the consensus. upstream_constant must begin at the
  # start codon; the tail of downstream_constant ends at the terminal stop.
  upstream_constant: ATGCAGCTT
  downstream_constant: ...GGATCCCATCACCACCACCATCACTAA
```

If your FASTA already holds complete assembled constructs, use `derive: true`
instead and declare nothing:

```yaml
    flanks:
      derive: true
      anchor_length: 20
```

**What changes in the output:** `index.csv` gains `insert_identity`,
`insert_edit_distance` and the two flank edit distances, so "the error is in the
vector, not my design" becomes a column rather than a guess. Acceptance floors are
measured over the insert, and chimera detection is scoped to it.

---

## 3. Pooled colony PCR, whole vector

**You have:** several culture plates combined into one 96-well colony-PCR plate, so
a well coordinate no longer identifies a source plate. What resolves it is the
gene: a design only ever picked onto plate 3 came from plate 3.

Add the pooling layout — this is **experimental design**, supplied by you, never
inferred from the reads:

```yaml
compressed_pcr:
  enabled: true
  block_pattern: "^(Block_\\d+)_"
  layout_csv: my_layout.csv
```

`my_layout.csv` is one row per (culture plate, block):

```csv
plate_barcode,culture_plate,library,block
BC01,PlateA,designs,Block_1
BC01,PlateB,designs,Block_2
BC01,PlateB,designs,Block_3
```

Check it before running:

```bash
nanopore3 layout --config my_run.yaml
```

which prints the layout and cross-checks it against your reference names, naming
any block with no plate or plate with no block.

**What changes in the output:** clones are filed under
`<barcode>/<culture plate>/<well>/` and `index.csv` gains a populated
`culture_plate`. A well that cannot be attributed to one plate keeps the
`|`-joined candidates and stays at well level — unresolved is reported as
unresolved, never guessed.

The complete version of this scenario is
[`configs/runs/260608_full_length_short.yaml`](../configs/runs/260608_full_length_short.yaml).

---

## Reading the results

`consensus_by_plate/index.csv` is the file most analyses should start from. One row
per clone:

| column | use |
|---|---|
| `grade` | `perfect`, `screenable`, `mixed_damage`, … — see [export.py](../src/nanopore3/export.py) |
| `culture_plate`, `well_id` | where to go back to on the bench |
| `insert_identity`, `insert_edit_distance` | accuracy of the designed part alone |
| `designed_allele_fraction` | for a mixed clone: how much of the well still carries the design |
| `mixed_worst_effect` | `silent` / `missense` / `nonsense`, or empty outside the ORF |
| `file` | the FASTA, relative to the tree root |

A `mixed_damage` clone is a real clone carrying a pre-transformation lesion, not an
artefact: usable, with the caveat those last two columns spell out.

---

## Changing your mind afterwards

Tuning a threshold does not mean repeating the demultiplexing, and does not need
the FASTQ still attached:

```bash
nanopore3 rerun --config tuned.yaml --from-run runs/<run-id> --from 04_consensus
```

Stages before `--from` are carried into a new run; everything after is recomputed.
If your change reaches back into an inherited stage, it is refused by name rather
than silently built upon. See [Architecture](architecture.md#reruns).
