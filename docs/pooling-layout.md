# Pooling layout: telling Nanopore3 which culture plate a well came from

Two experiments need two different amounts of configuration, and most runs are the
first kind.

## 1. One culture plate per barcode (monoclonal plates)

Nothing to configure. Leave `compressed_pcr` out of the configuration entirely, or
set `enabled: false`. A well is identified by its plate barcode and its well
barcode, and that is the whole answer.

```yaml
compressed_pcr:
  enabled: false
```

## 2. Several culture plates pooled into one colony-PCR plate

When several culture plates are combined into one 96-well colony-PCR plate, the
well barcode no longer identifies a source plate: well A1 of the PCR plate holds
A1 of every culture plate that went into it. What resolves the ambiguity is the
**gene**: if a well's read matches a design that was only ever picked onto culture
plate 3, that read came from culture plate 3.

That inference needs one table, and **that table is experimental design.** It is
what was done at the bench, it is known before any read exists, and Nanopore3 will
never infer it from the data. Inferring it would use the reads to construct the
key that then interprets those same reads, and the result would look like a
success no matter what was really on the plates.

### Supplying it as a spreadsheet (recommended)

```yaml
compressed_pcr:
  enabled: true
  layout_csv: my_layout.csv
  block_pattern: "^(Block_\\d+)_"
```

`my_layout.csv` is one row per **(culture plate, block)** pair:

```csv
plate_barcode,culture_plate,library,block
RP05,SUMO_A_P1,sumo_ab,A_Block_1
RP05,SUMO_A_P6,sumo_ab,A_Block_6
RP05,SUMO_A_P6,sumo_ab,A_Block_8
RP08,LAB_P1,lab,Block_1
```

Read it as one sentence per row: *this block was picked onto this culture plate,
and this culture plate went into this colony-PCR barcode.*

| column | meaning |
|---|---|
| `plate_barcode` | the colony-PCR plate's barcode, as it appears in `plate_reference_map` |
| `culture_plate` | the source plate you want recovered, named however you like |
| `library` | which reference library the block belongs to (block numbering usually restarts per library) |
| `block` | the assembly block, matching what `block_pattern` extracts from your reference names |
| `clonality` | optional: `per_block` if the design put one clone per block, else `unspecified`. Describes the **barcode**, so it must be the same on every row for that barcode. |

Repetition is expected and correct:

- a block split across several plates gets **one row per plate**
- a plate holding several blocks gets **one row per block**
- a duplicate row is ignored rather than counted twice

The file is resolved relative to the configuration file. A byte-order mark and
trailing blank lines - both of which spreadsheet exports produce - are tolerated.

### Supplying it inline

The same information can be written directly in YAML with `pcr_plates`, `blocks`
and `clonality`. Giving both forms is an error rather than a precedence rule: two
sources of truth for one design is how a stale block map becomes invisible.

### Migrating an existing configuration

```
nanopore3 layout --config my_run.yaml --export my_layout.csv
```

writes the current inline layout out as a table, then tells you what to put in the
configuration to use it. Run without `--export` to print a summary and check the
layout against the references.

## What is checked, and when

Preflight cross-checks the layout against the actual reference names and fails
before the run starts, naming every mismatch:

- blocks present in the references with **no culture plate** — reads matching them
  could never be attributed
- blocks in the layout matching **no reference** — usually a typo or a stale row
- culture plates that hold blocks but were **not pooled** into any barcode
- culture plates pooled into a barcode that hold **no block**
- barcodes in the layout missing from `plate_reference_map`

For example:

```
error: compressed_pcr layout does not match the references:
  library 'sumo_ab': 1 block(s) in the references have no culture plate: A_Block_7
  library 'sumo_ab': 1 block(s) in the layout match no reference: A_Block_99
  1 culture plate(s) were pooled but hold no block: GHOST_PLATE
```

The pooling design is separately checked for solvability when the configuration
loads: a layout in which two culture plates in one barcode share every block
cannot be deconvolved, and says so immediately rather than producing ambiguous
reads later.

## Where the block name comes from

By default, from the reference identifier, via `block_pattern`. The pattern's
first capturing group is the block name; with no group, the whole match is used.

```yaml
block_pattern: "^((?:[AB]_)?Block_\\d+)_"   # A_Block_12_dTF090_... -> A_Block_12
```

Check what it extracts before running:

```
nanopore3 layout --config my_run.yaml
```

Alternatively set `fragments_csv` on the reference library to take block identity
from a design table instead of the name.

## What you get back

Every consensus gains a `culture_plate` column, and the graded output tree files
clones under `<barcode>/<culture plate>/<well>/` instead of `<barcode>/<well>/`.
Where a read cannot be attributed to a single plate the field holds the
`|`-joined candidates and the clone stays at well level - an unresolved well is
reported as unresolved, never guessed.
