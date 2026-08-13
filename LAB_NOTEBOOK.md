# Nanopore3 lab notebook

A dated, append-only record of what changed, why, what was learned, and what to
do next. Newest entry first. Every entry should be readable on its own by a
person or an agent who has never seen this repository.

Conventions:

- Add a new entry rather than editing an old one. Corrections go in a new entry
  that names the entry it corrects.
- Record measurements with the command that produced them, not just the number.
- Separate **scientific** changes (which can alter results) from **performance**
  changes (which must not). Say explicitly which kind a change is.
- Reference evidence by path so it can be re-checked.

---

## 2026-08-13 (second) — Graded consensus tree, and `heterogeneous` renamed

### Graded per-plate output

Every run now writes `stages/05_qc/consensus_by_plate/`, one FASTA per consensus
under `<plate barcode>/<well>/`, named `<plate>_<well>__<design>__<grade>`:

```
RP06/A01/RP06_A01__Block_1_dTF083_156_2__perfect.fasta
RP06/A01/RP06_A01__Block_9_dTF083_282_1__frameshift.fasta
```

It lives under `05_qc` rather than `04_consensus` because the grade needs QC
results and a promoted stage directory is immutable. Wells are zero-padded so
`A01` sorts before `A10`, and names stay self-describing if a file is moved.

Full-run grades: **3,378 perfect**, 120 screenable, 124 frameshift, 121
mixed_variants, 35 truncated, 17 premature_stop, 7 mismatched, 2,072 low_depth.
89% of evaluable consensuses are exact matches to their design.

### `heterogeneous` renamed to `mixed_variants`

The experiment owner asked whether `heterogeneous` just means the well is
polyclonal, and if so to rename it. **It does not**, and the check is worth
recording:

- **468 of 512 wells hold more than one design** — 91% are polyclonal. That is
  the normal case here and is already expressed by such wells simply containing
  several files.
- Only **121 of 5,874 consensuses** are flagged.

The flag is set per *design*, when the reads assigned to one design disagree with
each other beyond the 0.60 support threshold, leaving ambiguity codes (median 5
bases, max 97). It means that one design is not a single clean clone: two
variants of it, cross-assigned reads, or noise.

The instinct was not baseless — flagged consensuses sit in wells with a median of
23 designs against 15 for clean ones, so crowded wells do produce them more
often. But naming it `polyclonal` would have conflated the expected condition of
91% of wells with a 2% anomaly.

Renamed to **`mixed_variants`** in both the grade vocabulary and the consensus
stage's `status` column, so one concept does not carry two names.

### A bug the rename exposed

Regenerating the tree after the rename left **121 stale `*heterogeneous*` files**
beside the new `*mixed_variants*` ones. The exporter wrote into an existing
directory without clearing it, so any rerun after a grade change, a threshold
change, or a re-analysis would leave files that no longer correspond to any
result. **A stale FASTA is indistinguishable from a current one once written**,
and would be read as a real screening result.

Fixed by building into a `.partial` sibling and swapping it in, matching the
atomic-promotion pattern the stages already use. An existing directory is only
replaced when it carries this exporter's `index.csv`, so an unrelated directory
is never deleted — tested both ways.

### Lessons

**"Rename X to Y" deserves the same check as any other claim.** Renaming to
`polyclonal` would have been a one-line change producing a permanently misleading
label on 121 files. Two counts — 468 of 512 wells versus 121 of 5,874
consensuses — settled it in a minute.

**A regenerable output must be regenerated destructively.** Writing over a tree
without clearing it is safe only while nothing is ever renamed or removed, which
is exactly the assumption that fails during a re-analysis.

### Next steps

Unchanged from the full-run entry, plus: consider surfacing the grade counts in
the HTML report and as a fourth figure, since they are now the most directly
useful summary of a run.

---

## 2026-08-13 — Full 2.17 M-read production run

**The full dataset has now been processed end to end.** Run
`runs/20260506-FULL-r1`, exit 0, 56.3 minutes wall clock, 2.0 GB of outputs.

### Stage timings

| Stage | wall time | throughput |
| --- | --- | --- |
| 01_ingest | <1 min | — |
| 02_demux | **45.9 min** | 786 reads/s |
| 03_assignment | **8.8 min** | 2,955 reads/s |
| 04_consensus | 1.3 min | — |
| 05_qc, 06_report | <1 min | — |
| **total** | **56.3 min** | |

Against an estimated 10+ hours before this work, but **demultiplexing came in
2.3x slower than the 1,818 reads/s benchmark** while assignment came in *faster*
than its 2,127 reads/s benchmark (startup amortises over 1.5 M reads). The
asymmetry is the explanation: assignment streams a local gzip, while
demultiplexing reads the 5.9 GB source FASTQ from **OneDrive cloud storage**. The
benchmarks were measured on a local file on an idle machine and were therefore
optimistic for this setting. Benchmark numbers should state where the input lived.

### Demultiplexing

| Outcome | reads | share |
| --- | --- | --- |
| assigned | 1,560,452 | 72.0% |
| well_unassigned | 243,063 | 11.2% |
| out_of_length | 191,675 | 8.8% |
| plate_unassigned | 170,590 | 7.9% |
| plate_ambiguous | 1,448 | 0.1% |
| well_ambiguous | 238 | 0.0% |
| **empty_read** | **91** | — |
| low_quality | 1 | 0.0% |

**`empty_read` is exactly 91**, matching the independent `awk` count made before
the feature existed. The reader accepts precisely the zero-length basecalls and
nothing else.

`low_quality` caught **one read in 2.17 million**, so the `minimum_mean_quality:
10` floor is doing essentially nothing on this data. Worth knowing before relying
on it as a control.

### Assignment, consensus, QC

| Assignment | reads | | Consensus (5,874) | | QC (5,874) | |
| --- | --- | --- | --- | --- | --- | --- |
| assigned_unique | 81.8% | | consensus_pass | 62.7% | pass | 60.5% |
| no_match | 12.9% | | low_depth | 35.3% | not_evaluable | 35.3% |
| motif_missing | 4.3% | | heterogeneous | 2.1% | fail | 4.2% |
| ambiguous | 1.0% | | | | | |

**93.5% of evaluable consensuses pass QC** (3,556 of 3,802), up from 88% on the
pilot, as expected once depth improves. `low_depth` remains 35% even at full
depth, so unlike the pilot artefact this is a real distribution: a third of
well/gene groups genuinely have fewer than six reads.

### Deconvolution validated at scale

**Zero `unexpected_block` across all 1,560,452 demultiplexed reads.** Every
assigned read's block is consistent with the declared layout. The failure mode is
loud — the earlier fictional layout produced 4,624 — so zero is strong evidence
the mapping is right.

RP01-RP04 report `unknown_pcr_plate` as intended; their layout is undeclared.

### Results worth acting on

**1. RP05 hits its prediction exactly.** Expected 6 clones per well from six
monoclonal culture plates mixed at colony PCR; **observed median is 6**.

**2. RP05 culture plate B1 has failed.** Share of resolved reads:

| | B1 | B2 | B3 | B4 | B5 | B6 |
| --- | --- | --- | --- | --- | --- | --- |
| RP05 | **0.8%** | 19.4% | 20.8% | 15.5% | 19.3% | 24.2% |

Five plates sit at 15-24%; block 1's plate contributes 0.8%. The pilot showed 1%
on 0.92% of the reads and the full run confirms it at depth. This looks like a
transformation or growth failure rather than a sequencing artefact.

**3. The two picking methods are distinguishable, as predicted.**

| | B1-3 | B4-6 | B7-10 |
| --- | --- | --- | --- |
| RP07 (monoclonal picks, mixed) | 37.2% | 31.8% | 30.9% |
| RP06 (scraped with a multichannel) | 41.2% | **22.7%** | 36.1% |

RP07 is near-even, as a deliberate equal mix should be. RP06 under-represents
blocks 4-6 by roughly a third. Scraping is measurably less uniform than mixing
monoclonal cultures.

**4. RP07 carries more diversity than designed.** Expected 10 clones per well
(3 + 3 + 4); **observed median is 14**. The pilot's 9 was an undercount from
sampling 0.92% of reads, and at full depth the wells exceed the design. So the
"monoclonal" picks feeding RP07 were not all monoclonal, or colonies carried
over during mixing. This is the clearest actionable discrepancy in the run.

### Lessons

**A prediction that lands exactly is worth more than a plausible number.** RP05's
median of 6 against an expected 6 validates the whole chain at once: barcode
calling, reference assignment, block attribution, and the layout model. It was
only possible because `clonality: per_block` derives an expectation from the
design instead of describing the data.

**Benchmarks inherit their storage.** Demultiplexing missed its benchmark by
2.3x purely because production input lives on cloud-synced storage. A throughput
figure without its input location is not reproducible.

### Next steps

1. **Check RP05 culture plate B1 at the bench.** 0.8% recovery is a failure, not
   a sampling artefact.
2. **Reconcile RP07's median of 14 against the designed 10.** Either the picks
   were not monoclonal or there was carryover during mixing.
3. Declare the RP01-RP04 layout if their culture-plate provenance matters;
   everything else about them is already processed.
4. Consider dropping `minimum_mean_quality` or lowering it, since it excluded one
   read in 2.17 million and is therefore not acting as a control.
5. Re-run the assembly-error detector (`scripts/detect_assembly_errors.py`) over
   the full demultiplexed reads to get a library-wide chimera and
   missing-fragment rate, now that a full run exists.
6. `low_depth` at 35% is real, not a pilot artefact. Decide whether
   `minimum_depth: 6` is the right threshold for the wells that matter.

---

## 2026-08-12 (ninth) — The real plate layout, derived rather than requested

### What changed

The previous entry ended by asking the experiment owner for the culture-plate
layout. That request was largely unnecessary: most of it had already been stated,
and the rest was derivable from the run itself. The layout is now committed in
both 20260506 profiles.

### What the data says about RP01-RP05

The stated model was "one block picked to one culture plate". The pilot
contradicts a literal reading of that. Blocks observed per plate barcode among
assigned reads:

| Plate | library | distinct blocks observed |
| --- | --- | --- |
| RP01 | aaseq_biotin | 18 of 20 |
| RP02 | aaseq_biotin | 20 of 20 |
| RP03 | aaseq_biotin | 18 of 20 |
| RP04 | aaseq_biotin | 19 of 20 |
| RP05 | dtf141_dtf142 | 6 of 6 |

Each of RP01-RP05 carries essentially its **whole library**, not a single block.
So these plates are not compressed at all: one culture plate maps to one colony
PCR plate, and gene identity has no plate ambiguity left to resolve. Their wells
are simply polyclonal.

**Compressed-PCR deconvolution therefore applies only to RP06 and RP07**, and for
those the layout was fully specified: three culture plates each, holding
`sumo_lab` blocks 1-3, 4-6 and 7-10.

### The committed layout

```yaml
compressed_pcr:
  pcr_plates:
    RP01: [CP_RP01]           # ... RP05 likewise, one culture plate each
    RP06: [CP_RP06_B1_3, CP_RP06_B4_6, CP_RP06_B7_10]
    RP07: [CP_RP07_B1_3, CP_RP07_B4_6, CP_RP07_B7_10]
  clonality:
    RP06: unspecified         # scraped: diversity per well is inconsistent
    RP07: per_block           # 3 + 3 + 4 = 10 expected clones per well
```

`sumo_lab` block 1 maps to **both** `CP_RP06_B1_3` and `CP_RP07_B1_3`. That is the
split-block case the validator was built for: it resolves because those two
plates sit in different colony PCR plates, so the reverse barcode separates them.

Culture-plate names are descriptive placeholders derived from the design. They
are identifiers, not measurements; rename them freely.

### Result

**Zero `unexpected_block` across all 14,984 reads.** Every assigned read's block
is consistent with the declared layout, which is the strongest available check
that the layout is right — a wrong mapping produces `unexpected_block` in bulk,
as the fictional demonstration layout did (4,624 reads). The only unresolved class
is `unknown_block` at 2,751, which is exactly the count of unassigned reads.

Source culture plate recovered, among resolved reads:

| Plate | B1-3 | B4-6 | B7-10 |
| --- | --- | --- | --- |
| RP07 (mixed monoclonal) | 36% | 32% | 32% |
| RP06 (scraped) | 41% | **23%** | 36% |

**This is a real result, not a QC number.** RP07 is near-even, as a deliberate
equal mix of monoclonal cultures should be. RP06 is skewed, under-representing
blocks 4-6 by a third — exactly what scraping a multichannel over an agar plate
would produce. The two picking methods are distinguishable in the output.

### Lessons

**Ask only for what cannot be derived.** The previous entry requested a layout
that was mostly already given and otherwise inferable from the block composition
of the reads. Deriving it first would have been faster and would have caught the
RP01-RP05 discrepancy sooner. Check what the data already knows before asking.

**A stated protocol and the observed data can disagree, and the disagreement is
information.** "One block per culture plate" does not describe RP01-RP05 as
sequenced. Rather than modelling around it silently, it is recorded here: either
the picking was not block-segregated, or the phrase meant something else. Either
way the analysis is unaffected, because those plates need no deconvolution.

**Zero is the informative number here.** `unexpected_block` at zero validates the
layout far more convincingly than any count of successes, because the failure mode
is loud: the fictional layout produced 4,624.

### Next steps

1. Run the full 2.17 M-read dataset; every blocker is now cleared and the layout
   is committed.
2. Confirm the RP01-RP05 reading above — is "one block per culture plate" a
   misremembering, or were those plates picked differently than described?
3. Re-measure the RP06/RP07 source-plate balance on the full run. The 23% dip for
   RP06 blocks 4-6 is from 0.92% of reads and needs the full depth before being
   quoted as a picking-efficiency result.
4. Compare observed clones per well against the expected 10 for RP07 once at full
   depth; the pilot's median of 9 is an undercount by construction.

---

## 2026-08-12 (eighth) — Polyclonal wells, per-library blocks, and figures

### Corrected understanding of the experiment

The experiment owner described the two picking modes, which changes the model:

- **RP01-RP05 (monoclonal picking, polyclonal PCR).** One block per culture
  plate, so deconvolution is trivial, but each well holds *many* genes from that
  one block. Diversity per well is real and unpredictable.
- **RP06 and RP07 (stacked polyclonal).** Three culture plates each, holding
  blocks 1-3, 4-6 and 7-10, pooled into one colony PCR plate. **RP06** was
  scraped with a multichannel over an agar plate, so its per-well diversity is
  inconsistent. **RP07** came from monoclonal picks mixed and grown together, so
  it should carry close to 3, 3 and 4 clones from the three source plates: ten.

The pilot data agrees before any expectation was configured. Distinct genes per
well, over 0.92% of reads so certainly an undercount:

| Plate | wells | median | max |
| --- | --- | --- | --- |
| RP01 | 54 | 3 | 6 |
| RP03 | 25 | 1 | 10 |
| RP05 | 94 | 5 | 6 |
| RP06 | 95 | **8** | 14 |
| RP07 | 95 | **9** | 13 |

RP07's median of 9 against an expected 10 is exactly the predicted shape, and
RP06 is both lower and more dispersed, as scraping implies.

### A bug this exposed

**Block identifiers collide between reference libraries.** `aaseq_biotin` has
blocks 1-20, `dtf141_dtf142` 1-6, `sumo_lab` 1-10. The map committed in the
previous entry was keyed on the block alone, so block "1" of one library would
have silently resolved using another library's culture plate. `blocks` is now
nested by library and `resolve()` takes the library. Caught only because this
experiment routes three libraries; a single-library run would never have shown it.

### What was added

- **`clonality` per colony PCR plate.** `per_block` derives the expected clones
  per well from the layout — the blocks across that plate's pooled culture
  plates. `unspecified` (the default) reports observed counts only. A scraped
  plate has no expectation, and inventing one would manufacture deviations.
- **`src/nanopore3/figures.py`** and `scripts/make_figures.py`. Every run now
  writes figures into `stages/06_report/figures/`, best-effort so a missing
  matplotlib never fails a run; `figures.json` records what happened.

### Figure design decisions

Sized for print: 89 mm single / 183 mm double column, 6-8 pt sans, vector PDF plus
600 dpi PNG, `pdf.fonttype 42` so text stays editable.

**The palette is computed, not chosen.** `scripts/validate_palette.py` implements
the checks in Python — OKLab conversion, Machado-Oliveira-Fernandes (2009)
protanopia and deuteranopia at severity 1.0, pairwise dE, lightness band, chroma
floor and contrast — because the environment has no Node to run the reference
validator, and colour-vision safety should not be judged by eye.

Measured results, which decided the palette:

- `#0072B2, #D55E00, #009E73` passes all-pairs: worst CVD dE **11.0**, worst
  normal-vision dE 18.7.
- Adding `#CC79A7` drops the worst CVD pair to **7.6**, below the target of 8, so
  three is the cap wherever two marks can touch.
- The common Okabe-Ito orange `#E69F00` **fails contrast** at 2.25:1 on white,
  under the 3:1 floor for marks; `#D55E00` replaces it at 3.9:1.

Layout problems found only by rendering and looking, which is why that step is
not optional:

1. Plate-map tick labels collided with the panel titles below them. Fixed by
   drawing coordinates only on the outer edge — they repeat in every panel, so
   96 repetitions per panel is noise.
2. Two columns of 96-well panels forced a figure taller than a page. The panel
   aspect is fixed at 12:8, so the figure height is now derived from the column
   count rather than guessed, and the spare grid cell holds the colour bar
   instead of being left as a hole.
3. The clonality legend sat on top of RP07's data, and advertised an "expected"
   series that was not drawn because no plate declared one. The legend now sits
   outside the axes and only lists series actually present.

### Lessons

**A layout bug is invisible to a validator.** The palette check passed on the
first attempt; every real defect was geometric and only appeared on screen.
Render and look, every time.

**Defaults must not invent expectations.** `clonality` defaults to
`unspecified`, so a plate only gets an expected-diversity line when the design
actually justifies one. A default of "one clone per block" would have drawn a
confident reference line on a scraped plate and turned an unknown into an
apparent failure.

**A collision that only appears with three libraries.** The block-key bug was
latent and would have produced confident wrong culture plates rather than an
error. Multi-tenant keys need qualifying the moment a second tenant exists.

### Next steps

1. **Supply the real culture-plate names and block groupings** for RP01-RP07.
   Everything here is built and tested but has still only run against a
   demonstration layout.
2. Set `clonality: per_block` for RP07 and `unspecified` for RP06 once those
   names exist, then compare observed against expected per well.
3. RP01-RP05 diversity is currently undercounted because the pilot samples 0.92%
   of reads; re-measure on the full run before drawing conclusions about picking
   efficiency.
4. Consider a per-well figure of *which* source culture plate was recovered, to
   expose systematic dropout of one plate in the stack.

---

## 2026-08-12 (seventh) — Compressed-PCR deconvolution

### The workflow being modelled

Colonies from one assembly block are picked into a culture 96-well plate. Several
culture plates are then pooled into **one** colony-PCR plate, which carries the
forward (well) and reverse (plate) primer barcodes. So a barcode pair identifies a
*PCR* well, which may hold colonies from several culture plates at the same
position.

The assigned gene supplies the missing coordinate. A gene belongs to exactly one
block, and a block was picked into known culture plates, so
``gene -> block -> culture plate`` recovers the source. When a block is split
across culture plates, those plates must go to **different** colony PCR plates so
the reverse barcode separates them.

### What was built (scientific, opt-in)

`src/nanopore3/deconvolution.py` with a `compressed_pcr` config section:

```yaml
compressed_pcr:
  enabled: true
  pcr_plates:            # culture plates pooled into each colony PCR plate
    RP01: [CP_A, CP_B, CP_C]
    RP02: [CP_D]
  blocks:                # culture plate(s) each block was picked into
    "3": [CP_C, CP_D]    # split, but across different PCR plates
```

Design decisions worth recording:

- **An unresolvable layout is rejected at configuration time.** If two culture
  plates holding the same block are pooled into one PCR plate, no evidence can
  separate them. That is detectable statically, so it is a config error with a
  message naming the offending plates, rather than ambiguous reads discovered
  after a multi-hour run.
- **Ambiguity is never resolved arbitrarily.** If several plates remain
  candidates the read is reported `ambiguous` with no plate assigned.
- **`unexpected_block` is a real signal, not an error state.** A gene whose block
  was never pooled into that PCR plate indicates a mis-pick, cross-contamination,
  or a misdescribed layout. It is worth watching.
- **The consensus grouping key is deliberately unchanged.** `culture_plate` is
  recorded alongside consensus rows and FASTA headers as provenance, so consensus
  identities stay byte-identical whether or not deconvolution runs. Folding it
  into identity would have silently changed every prior consensus ID.

### Where the block comes from

Two independent sources, which cross-validate:

- the reference identifier, via `block_pattern` (default `^Block_(\d+)_`);
- the oPool design table, via `fragments_csv` per library, joined **by sequence**
  because the FASTA and design table use different identifiers.

Checked across all three libraries: **1,175 of 1,175 references agree between the
two sources, zero disagreements.** The design table is treated as authoritative
when supplied.

### Verification

A demonstration layout over the 20k pilot (7 PCR plates, 20 culture plates)
resolved reads to 5-6 distinct source culture plates per PCR plate, which is the
deconvolution working. Status counts were `resolved` 7,609, `unexpected_block`
4,624, `unknown_block` 2,751 out of 14,984.

The two non-resolved classes are exactly right for that input and worth
explaining, because both look alarming:

- `unknown_block` 2,751 is **exactly** the number of unassigned reads
  (14,984 - 12,233 = 2,751). No gene means no block; nothing is wrong.
- `unexpected_block` 4,624 is an artefact of the **invented** demonstration
  layout, which mapped block *N* to culture plate *N* in every library. Real
  reads therefore landed in PCR plates that layout says never held their block.
  This is the mis-layout detector working, on a deliberately wrong layout.

**No real layout has been supplied yet**, so the demonstration config was not
committed; a plausible-looking but fictional plate map is exactly the kind of
artefact that later gets mistaken for real metadata.

### Also fixed

`rescue_policy` and `rescue_candidates` existed on `ReferenceSettings` and were
documented as per-library configuration, but were never added to the YAML parser,
so a run file setting them would have been **rejected as an unknown key**. They
are now parseable. Documentation claiming a knob exists is worse than no knob.

### Next steps

1. **Supply the real block-to-culture-plate map** and re-run. Until then this
   feature is untested against a true layout.
2. Watch `unexpected_block` on the real layout: a non-trivial rate means
   mis-picking or cross-contamination between culture plates.
3. Consider reporting resolved reads per culture plate well in the HTML report,
   so plate-level dropout is visible at a glance.
4. Deconvolution currently keys on the block alone. If a future design puts two
   blocks in one culture plate, the mapping still works; if one block spans more
   culture plates than there are PCR plates, it cannot, and the config error will
   say so.

---

## 2026-08-12 (sixth) — CORRECTION: blocks bound recombination, not overhangs

**This entry corrects the 2026-08-12 (fifth) entry. Every chimera call it
reported was a false positive.**

### What was wrong

The fifth entry grouped interchangeable fragments by overhang pair across a whole
reference library. The experiment owner corrected the model: `Block` is the
**sub-pool** a gene is assembled in. Golden Gate overhangs are unique *within* a
block and deliberately **reused between** blocks, which is safe because each block
is separately PCR-amplified with its own primer pair and ligated in its own
reaction. Two blocks never share a tube.

Confirmed in the design tables — every block has exactly **one** primer pair:

| Library | genes | blocks | within-block overhang pairs | shared within a block |
| --- | --- | --- | --- | --- |
| aaseq_biotin | 651 | 20 | 1,041 | 2 |
| dtf141_dtf142 | 192 | 6 | 259 | 1 |
| sumo_lab | 332 | 10 | 545 | 1 |

Within a block, junction overhangs are **unique**: 1,039 of 1,041 pairs have a
single variant. The handful of exceptions are the vector pair `GCTT`/`AGTG`
shared by single-fragment genes, which is not a junction at all.

So the previous model's "4–18 interchangeable variants per overhang pair" were
genes from *different blocks*. All 15 reported chimeras were checked against the
block assignments: **0 same-block, 15 cross-block — every one physically
impossible.** The number was not merely imprecise, it was meaningless.

### The corrected model

Two consequences, and the second is the interesting one:

1. `FragmentLibrary` is now block-aware. `interchangeable()` requires a block and
   `Gene` carries its block. Overhang identity alone no longer implies anything.
2. **A chimera cannot be a legitimate overhang swap.** Since within-block
   overhangs are unique, correct ligation cannot mis-pair. Real chimeras must come
   from mis-ligation of *similar but non-identical* overhangs, or from PCR
   template switching during block amplification. Either way the recombination
   space is **the whole block**, not the overhang class.

The detector therefore now compares each observed fragment against the
corresponding slot of **every gene in the same block**.

### Corrected result on the same 3,000 pilot reads

| Outcome | before (wrong) | after (block-scoped) |
| --- | --- | --- |
| intact | 1,819 | 1,819 |
| single-fragment gene | 581 | 581 |
| unresolved | 241 | 232 |
| fragment_missing | 181 | 181 |
| no insert | 163 | 163 |
| **chimeric** | 15 (all bogus) | **24** |

**24 chimeric reads, 1.3% of decided reads, and all 24 have same-block donors.**
Per-slot identities are 0.91–0.98, so each fragment confidently matches a
different gene in the same reaction. Examples:

```
block  3  dLK27_325_c3   >> binder_dLK38_422_model   identity 0.927 / 0.978
block  3  dTF142_19_c5   >> dTF141_91_c5             identity 0.963 / 0.975
block  7  dTF082_dTF080_l125_... >> dTF083_152_4     identity 0.926 / 0.912
block  5  dTF083_451_3   >> dTF086_432_1
```

The last two independently reproduce pairs found by the crude split-half probe in
the fourth entry, which reported `Block_5 dTF083 + dTF086` and
`Block_7 dTF082 + dTF083`. Two methods with different assumptions agreeing on the
same recombinants is the strongest evidence so far that these are real.

### Lessons

**A biological constraint I did not know about silently defined the answer.** The
overhang model was internally consistent, reconstructed all 1,175 genes exactly,
and produced confident, plausible, named chimeras. It was still entirely wrong,
because physical separation into sub-pools — a fact recorded in a column I had
read but treated as a label — determines what can recombine. No amount of
internal validation would have caught it; only the domain owner could.

**"Plausible and named" is not "verified".** The false calls looked more credible
than the true ones because they came with donor names and identities. The check
that mattered was trivial once the constraint was known: are the two donors in the
same block? It took one query and invalidated the whole result.

**Ask what makes a hypothesis physically impossible, not just unlikely.** The
useful question was not "how similar are these sequences" but "were these two
molecules ever in the same tube".

### Next steps

1. Review a few of the 24 same-block chimeras by eye before trusting the class.
2. Check whether chimeric donors have *similar* junction overhangs, which would
   distinguish mis-ligation from PCR template switching. If overhangs differ
   greatly, template switching is the likelier mechanism.
3. Consider reporting per-block assembly quality: blocks with more recombination
   may have poorly separated overhangs.
4. Unchanged from before: 19.4% of reads are single-fragment genes and cannot be
   assessed for mis-assembly by this method at all.

---

## 2026-08-12 (fifth) — Golden Gate fragment model and assembly-error detection

### What changed

The experiment owner supplied the oPool design tables and confirmed the critical
detail: **the oligo sequences in the pool are not what ends up in the gene.** Each
oligo carries a BsaI site and primer padding at both ends.

**1. New module `src/nanopore3/fragments.py`.** Models the Golden Gate assembly:

- Excises each oligo at `GGTCTC N <4 nt overhang> ... <4 nt overhang> N GAGACC`,
  keeping the overhangs and discarding the recognition sites and padding.
- Joins neighbours *through* their shared overhang, so a junction overhang
  appears once in the product, and strips the two outermost overhangs, which are
  contributed by the vector.
- Validates every gene by reconstruction and **raises rather than guessing**. A
  silently wrong fragment model would misreport every assembly call downstream.
- Groups fragments by overhang pair. Two designs can only mis-assemble where they
  share a pair, so those variants are exactly the set a chimera is drawn from.

A subtlety worth recording: **taking the first `GGTCTC` in an oligo is wrong.**
Primer padding sometimes contains one, which silently produced an over-long
fragment. Reconstruction rates were 174/192 and 642/651 with naive selection.
Choosing the BsaI pair that actually reconstructs the gene gives **100% on all
four design tables**.

**2. `scripts/detect_assembly_errors.py`.** Projects a gene's designed fragment
boundaries onto each read through an alignment, then re-matches each observed
fragment against the interchangeable variants for that slot's overhang pair.

### Which SUMO design table is correct

The experiment owner was unsure between two candidates. Resolved by matching
sequences against the `sumo_lab` reference FASTA:

| Candidate | sequences matching the 332 references |
| --- | --- |
| `opTF001/dTF017_..._FULL_INFO.csv` | **332 / 332** |
| `opTF005/LAB/LAB_FULL_INFO.csv` | 26 / 332 |

**`opTF001/dTF017_dTF018_dTF020_dTF083_dTF084_dTF086_dTF082_dTF023_dTF024` is the
correct table.** Note that all three libraries match by *sequence* but **0 by
name** — the reference FASTA identifiers carry a `Block_N_` prefix that the design
tables do not. Joins between the two must be sequence-based.

### Design structure

| Library | genes | fragments/gene | overhang pairs | interchangeable |
| --- | --- | --- | --- | --- |
| aaseq_biotin | 651 | 1:135, 2:509, 3:7 | 72 | 65 |
| dtf141_dtf142 | 192 | 1:63, 2:129 | 65 | 65 |
| sumo_lab | 332 | 1:60, 2:272 | 65 | 65 |

Fragments are ~205 nt median. Interchangeable sets hold 4-18 variants, plus one
large class per library (60-135) formed by the single-fragment genes, which all
share the vector overhang pair `GCTT`/`AGTG`.

`VectorOH1` is `GCTT` and `VectorOH2` is `AGTG`, independently confirming the
boundary work in the earlier entries: `ATGCAGCTT` ends with the 5' vector
overhang and `AGTGGATCC` begins with the 3' one.

### Result on 3,000 pilot reads

| Outcome | reads | share |
| --- | --- | --- |
| intact | 1,819 | 60.6% |
| single-fragment gene (no junction to mis-pair) | 581 | 19.4% |
| unresolved | 241 | 8.0% |
| fragment_missing | 181 | 6.0% |
| no insert extracted | 163 | 5.4% |
| **chimeric** | 15 | 0.5% |

**0.8% of reads with a decided outcome are chimeric**, with donors named, e.g.
`dTF142_199_c5 >> dTF141_23_c5` and `dLK10_719_1 >> dLK10_862_1` — recombination
between designs, across the two dTF14x families in the first case.

This supersedes the crude split-half figure in the previous entry. That probe
reported 48% discordance *among resolvable `no_match` reads only*, which is a
small and heavily selected subset; it is not a library-wide chimera rate. The
fragment-aware number is the one to quote.

`fragment_missing` at 6% is a real and separate assembly failure class.

### Lessons

**The reagent is not the product.** The single most important fact here was that
an oligo contains a BsaI site and padding that never reach the assembled gene.
Any analysis matching reads against raw pool oligos would be wrong by ~40 nt per
fragment, and would have looked plausible.

**Validate a derived model by reconstruction, not by inspection.** The naive BsaI
excision looked correct on the first row examined and was wrong on 9% of genes.
Requiring every gene to reconstruct exactly turned a plausible parser into a
verified one and located the padding-site problem immediately.

### Next steps

1. **Decide whether to promote this into the pipeline** as a `chimera` state in
   `03_assignment`, which `docs/architecture.md` already reserves. It is
   deliberately a separate script for now, so the validated pipeline outputs do
   not change before the calls have been reviewed.
2. **Review a handful of the 15 chimeric reads by eye** before trusting the class.
3. Add fragment CSV paths to the run profile once integrated, keyed per library.
4. Investigate the 8% `unresolved` and 6% `fragment_missing`; some will be
   nanopore noise on ~205 nt fragments rather than genuine assembly failures.
5. Chimera detection currently needs >=2 fragments, so 19.4% of reads (single
   fragment genes) cannot be assessed this way at all. For those, mis-assembly
   would have to be detected as recombination *within* a fragment, which the
   overhang model cannot see.

---

## 2026-08-12 (fourth) — Assembly failure modes: chimeras are not detected

**No code changed in this entry.** It records a measured gap, so the next person
does not assume assignment handles oligo-pool assembly errors.

### Context

The experiment owner clarified that this is Nanopore sequencing of an **oligo-pool
assembly**. Molecules are physically assembled from fragments, so the expected
failure modes are a **missing fragment**, an **extra fragment**, and
**mis-pairing** — a chimera whose fragments come from different references.

### How the current pipeline actually behaves

Constructed from two real 372 nt `aaseq_biotin` references and run through the
production thresholds (identity 0.80, coverage 0.70/0.70, margin 0.02):

| Construct | Assignment | QC (assuming a perfect consensus) |
| --- | --- | --- |
| intact | `assigned_unique`, identity 1.000 | pass |
| 1/6 fragment missing | `assigned_unique`, identity 0.833 | **fail** (length, frame) |
| 1/3 missing | `no_match`, ref coverage 0.637 | — |
| 1/2 missing | `no_match`, ref coverage 0.468 | — |
| 1/6 extra inserted | `assigned_unique`, identity 0.857 | **fail** (length, frame) |
| 1/3 extra inserted | `assigned_unique`, query coverage 0.742 | **fail** (length, frame) |
| chimera 50% A + 50% B | `no_match`, margin 0.005 | — |
| chimera 67/33 | **`assigned_unique` → A**, margin 0.175 | fail (identity only) |
| chimera 75/25 | **`assigned_unique` → A**, margin 0.250 | fail (identity only) |
| chimera 90/10 | **`assigned_unique` → A**, identity 0.941, margin 0.392 | fail (identity only) |

**Indels are handled acceptably; chimeras are not.** A missing or extra fragment
changes the length, so `expected_length` and `reading_frame` catch it
independently of identity. A chimera of two same-length designs is the *right
length and in frame*, so those checks pass. The only thing standing between a
chimera and a clean report is the 0.98 identity gate — the same scalar that
low-depth consensus noise moves. A 90/10 chimera scores 0.941, which is not
cleanly separable from a noisy but correct consensus.

Worse, the margin check gives no protection. It compares the best against the
second best *whole-length* alignment, so an uneven chimera looks *more* confident
the more lopsided it is: margin rises from 0.175 at 67/33 to 0.392 at 90/10.

### How much of this is in the real data

A split-half probe over 4,000 pilot reads: cut each extracted insert in half,
shortlist each half independently by specificity-weighted k-mers, and compare.

| Assignment status | halves agree | **halves differ** | unresolved |
| --- | --- | --- | --- |
| `assigned_unique` | 3,214 | **24 (0.7%)** | 1 |
| `no_match` | 90 | **83 (48%)** | 269 |
| `ambiguous` | 0 | **3 (100%)** | 0 |

- **Roughly half of all resolvable `no_match` reads are chimeric.** They are
  currently discarded into the same bucket as genuine junk, so a real and
  measurable assembly failure rate is being thrown away rather than reported.
- **0.7% of `assigned_unique` reads are chimeric and reported as clean.** Small
  as a fraction, but it is a *confident wrong answer*, and it scales to thousands
  of reads across 2.17 M.

The discordant pairs are highly informative:

```
5' Block_5_dTF083_451_3      3' Block_5_dTF086_432_1
5' Block_7_dTF082_dTF080_l1  3' Block_7_dTF083_152_4
5' Block_9_dTF017_6i2g_l91_  3' Block_9_dTF024_APdesign_
```

**The two parents share the same `Block_N` prefix every time.** Mis-pairing
happens between designs occupying the same assembly block, which is exactly what
a shared fragment junction predicts. The reference naming already encodes the
block structure needed to model this.

### Recommended approach

Detection does not need new algorithms — the existing specificity-weighted k-mer
index already resolves each half independently, as the probe above shows. What is
needed is to run it **segment-wise** and record the profile:

1. Split the extracted insert into N windows (or, better, at the known fragment
   boundaries from the oPool design).
2. Shortlist each window against the routed library.
3. Uniform best reference across windows and a passing whole-length alignment →
   `assigned_unique`, as now. Windows confidently naming different references →
   a new `chimera` state carrying both parents and the breakpoint window.
4. Require a margin per window so noisy windows abstain rather than vote.

`docs/architecture.md` already reserves `chimera` and `truncated` in the
assignment vocabulary, and the README already states that v0.3 does not make
biological chimera calls — so this closes a documented gap rather than inventing
a contract.

Using the real fragment boundaries would be materially better than fixed windows,
because it would also distinguish *missing* and *extra* fragments by name rather
than inferring them from length. That needs the oPool design (fragment sequences
or coordinates per reference), which is not in this repository.

### Next steps

1. Decide between fixed-window and fragment-aware segmentation; the latter needs
   the oPool fragment definitions.
2. Implement segment-wise assignment and the `chimera` state.
3. Report an assembly-failure breakdown per plate/well: correct, missing
   fragment, extra fragment, mis-paired. This is a scientific result about the
   assembly, not merely a QC filter.
4. Until then, **treat `no_match` as "rejected, cause unknown" and do not quote
   it as a contamination or quality figure** — about half of it is chimeric.

---

## 2026-08-12 (third) — Coding QC, and the 6-base offset closed

### Context

The experiment owner confirmed the open question from the previous entry: the
first nine bases of every construct are **`ATGCAGCTT`**, so `CAGCTT` is a constant
linker between the start codon and the variable insert. They also supplied the
constant region 3' of the insert (555 nt, beginning `AGTGGATCC`, ending in a
His-tag and `TAA`), and asked for reading-frame and internal-stop QC that adapts
to whatever constant regions a user declares.

### What changed

**1. Reading-frame and internal-stop QC now actually run (scientific).**
Previously both were permanently `not_evaluable`: they required
`coding_start`/`coding_end` integers that no configuration could supply.

They are now driven by two declared constant regions in the `qc` section:

```yaml
qc:
  upstream_constant: ATGCAGCTT     # must begin at the start codon
  downstream_constant: AGTGGATCC...TAA
```

QC assembles `upstream_constant + consensus + downstream_constant`, requires the
total to be a whole number of codons, translates frame 0, and fails if a stop
appears before the final codon. Design decisions worth knowing:

- **Both constants are required together.** A frame inferred from one side would
  be a guess, and a stop-codon result computed in a guessed frame is worse than
  no result. Supplying one without the other is a configuration error.
- **`upstream_constant` must begin at a start codon** (`ATG`/`GTG`/`TTG`),
  validated at load time. This is what makes the frame declared rather than assumed.
- **Stops are not reported from a frameshifted ORF.** When the frame check fails,
  `internal_stops` is `not_evaluable`, because stop codons read out of frame are
  noise and would be reported as if they were biology.
- **Ambiguous codons translate to `X`, never to a guessed residue**, so a
  low-support consensus base cannot become a confident amino acid.
- Omitting both keeps the previous behaviour exactly, so this is opt-in.

Two new columns in `qc.csv.gz`: `protein_length` and `internal_stop_codon`.

**2. The 6-base offset is closed (scientific).**
With `CAGCTT` confirmed constant, `forward_motif` is extended through it to
`TAAGAAGGAGAGCAGCTATGCAGCTT`, and `motif_max_edits` raised 2 → 3 because a 26 nt
motif at 2 edits is proportionally stricter than the original 20 nt motif was.
`qc.upstream_constant` is correspondingly `ATGCAGCTT`, since the consensus no
longer carries the linker.

**The median consensus-to-reference edit distance is now 0, previously 6.**

### Results on the 20k pilot

| | at session start | + length gate | + coding QC | + motif fix (final) |
| --- | --- | --- | --- | --- |
| `assigned_unique` | 12,145 | 12,141 | 12,141 | **12,233** |
| `motif_missing` | 973 | 899 | 899 | **668** |
| `ambiguous` | 107 | 50 | 50 | 167 |
| QC pass | 286 | 286 | 273 | **663** |
| QC fail | 466 | 466 | 479 | 91 |

Final QC breakdown over 754 evaluable consensuses:

| Criterion | pass | fail | not_evaluable |
| --- | --- | --- | --- |
| `full_amplicon` | 680 | 74 | — |
| `expected_length` | 737 | 17 | — |
| `reading_frame` | 730 | 24 | — |
| `internal_stops` | 726 | 4 | 24 |

Protein length: min 232, median 280, max 347 residues.

**QC passes went from 286 to 663 while two additional real checks were added.**
On the pre-motif-fix run the coding checks alone caught **13 consensuses that
passed both identity and length** but were frameshifted or carried a premature
stop — defects the previous QC could not have reported.

### Lessons

**A systematic offset masquerades as a quality problem.** Every symptom pointed at
consensus quality: identity just under threshold, a plausible-looking ~2%
error rate, a pass rate that looked depth-limited. The actual cause was a constant
6-base boundary disagreement between motif extraction and the reference
definition. The tell was that median identity was *flat across every depth
bucket* — real sequencing error improves with depth, systematic offsets do not.
Check whether an error metric responds to the variable that should govern it.

**A quantised error distribution is not sequencing error.** Median edit distance
of exactly 6, with p10 also exactly 6, cannot come from a stochastic process. The
distribution's shape identified the bug before any alignment was inspected.

**Half a QC suite silently not running is worse than not having it.** The stage
reported `reading_frame` and `internal_stops` for every consensus, always as
`not_evaluable`, because nothing could supply their parameters. It looked like
coverage. For a protein-design assay these are the checks that matter most, and
nothing was frameshift-checked until now.

### Next steps

1. **Run the full 2.17 M-read dataset.** All blockers are now cleared and the
   profile is settled. Compare the `empty_read` count against the expected 91.
2. **Review the 24 frameshifted and 4 premature-stop consensuses** — these are
   real biological or synthesis defects, and are the first such calls this
   pipeline has ever made. Confirm a few by eye before trusting the class.
3. **Investigate the remaining 74 `full_amplicon` failures** now that the
   systematic offset is gone; whatever remains should be genuine variation.
4. `ambiguous` rose 50 → 167 with the relaxed motif budget. These are reads
   rescued from `motif_missing` and conservatively flagged rather than newly
   confused, but confirm that reading before quoting yields.
5. Rename or fix the length-ratio "coverage" metrics in `qc.py` (carried over).

---

## 2026-08-12 (later) — Read-length gate, and the QC failures are a 6-base offset

### What changed

**1. Read-length gate of 500–2,000 nt (scientific, requested).**
The experiment owner supplied a conservative range to exclude fragments and
primer dimers. Set in both 20260506 profiles as `library.minimum_read_length`
and `maximum_read_length`.

Calibration on the 20k pilot: assigned reads span **974 nt (p1) to 1,462 nt
(p99)**, median 1,247, so 500–2,000 sits well outside the real distribution and
is genuinely conservative. It removes 1,021 of 20,000 reads (5.1%), but **86% of
those were already failing demultiplexing** — only 142 reads (0.94%) that would
otherwise have been assigned are affected.

Effect on the pilot:

| | no gate | 500–2,000 |
| --- | --- | --- |
| demux assigned | 15,126 | 14,984 |
| assignment `assigned_unique` | 12,145 | 12,141 |
| assignment `ambiguous` | 107 | **50** |
| assignment `motif_missing` | 973 | 899 |
| `consensus_pass` | 603 | 603 |
| QC pass | 286 | 286 |
| demux + assignment wall time | 30 s + 20 s | **19 s + 13 s** |

Four assigned reads lost, assignment ambiguity **halved**, and the run about a
third faster. The gate is close to free scientifically and clearly worth it.

**2. Read-level gates now short-circuit before barcode matching (performance).**
`_demux_one` computed `length_status` and `quality_status` up front but still ran
both barcode panels before applying them, so a length filter cost full compute.
Length and quality already override any barcode verdict, so a failing read now
returns immediately with `plate`/`well` recorded as `not_attempted`. This is what
makes the filter cheap; it is also why filtered reads no longer carry barcode
evidence columns. No previously-configured run had a length gate, so no existing
result changes.

### The QC failures are not a depth problem — they are a constant 6-base offset

The previous entry's next-steps list called the 1,731 `low_depth` groups and 466
QC failures a "reads-per-well problem". **Both halves of that were wrong.**

**`low_depth` is an artefact of the pilot's size.** The pilot is 20,000 of
2,167,558 reads (0.92%). Median depth is 3 reads per group and only 752 of 2,483
groups reach the `minimum_depth: 6` threshold. At full scale each group receives
roughly 100× more reads. `minimum_depth: 6` is unchanged and appropriate; nothing
needs adjusting.

**The QC failures are a systematic boundary offset.** QC pass rate does not
improve with depth — it plateaus around 48% and median identity is flat at
~0.978 in every depth bucket from 6 to 30+. Aligning each consensus to its
reference shows why:

```
cigar: 2I1=1I1=3I238=      cons 5' CAGCTTGTCGTTGGTGGAGTG...
                           ref  5' ------GTCGTTGGTGGAGTG...
```

Six inserted bases at the 5' end, then a **perfect** match for the entire
remainder. **707 of 752 consensuses (94%) carry the identical prefix `CAGCTT`
that the reference sequences do not include**, across all three libraries.

This also explains the pass/fail split, which is otherwise puzzling. Identity is
`1 - 6/length`, so a reference of 300 nt or more still clears the 0.98 gate with
those 6 edits while a shorter one cannot. **QC outcome was being determined by
insert length, not by consensus quality.** The `minimum_identity: 0.98` threshold
is not too strict; the sequences being compared were misaligned by 6 bases.

Two diagnostic runs testing `forward_motif: TAAGAAGGAGAGCAGCTATG` **+ `CAGCTT`**:

| | baseline | +CAGCTT, edits 2 | +CAGCTT, edits 3 |
| --- | --- | --- | --- |
| `assigned_unique` | 12,141 | 12,011 | **12,233** |
| `ambiguous` | 50 | 138 | 167 |
| `motif_missing` | 899 | 968 | **668** |
| QC pass | 286 | 665 | **680** |
| QC fail | 466 | 74 | **74** |

Extending the motif alone costs assignment yield, because a 26 nt motif at
`motif_max_edits: 2` is proportionally stricter than a 20 nt one. Raising the
budget to 3 recovers it and then some: against the baseline it rescues 231 reads
that were discarded as `motif_missing`, assigns 92 more, conservatively flags 117
more as `ambiguous`, and **more than doubles QC passes**.

**This change was deliberately not applied in this entry.** It depended on a fact
only the experiment owner could confirm: whether `CAGCTT` is a constant linker
between the ATG and every designed insert — in which case trimming it is correct
and lossless — or whether it varies for some designs, in which case trimming would
corrupt them. Evidence is preserved in `runs/20260506-diag-motif-cagctt` and
`runs/20260506-diag-motif-cagctt-e3`.

> **Resolved.** The experiment owner confirmed `ATGCAGCTT` as the constant first
> nine bases of every construct. Applied in the 2026-08-12 (third) entry above.

The alternative fix is to prepend `CAGCTT` to the reference FASTAs instead, which
keeps the linker in the consensus. Both make consensus and reference describe the
same molecule; extending the motif is a pure configuration change and is cheaper.

### What QC actually evaluates

`evaluate_consensus` in `src/nanopore3/qc.py` scores four independent criteria and
fails `overall` if any **evaluable** one fails:

| Criterion | Test | Status here |
| --- | --- | --- |
| `full_amplicon` | global NW identity ≥ 0.98, and both length-ratio coverages ≥ 0.95 | the only one that ever fails |
| `expected_length` | `abs(len(consensus) - len(reference)) <= 10` | passes 738/752 |
| `reading_frame` | coding length divisible by 3 | **never runs** |
| `internal_stops` | no TAA/TAG/TGA before the final codon | **never runs** |

Two gaps worth knowing about:

- **Half of the advertised QC never executes.** `reading_frame` and
  `internal_stops` require `coding_start`/`coding_end`, and the pipeline never
  passes them — there is no configuration field for them at all. They are
  permanently `not_evaluable`. For a protein-design assay these are the checks
  that matter most, so this is a real gap, not a cosmetic one.
- **`query_coverage` and `reference_coverage` are length ratios, not coverage.**
  `global_alignment_metrics` computes `min(1, len(ref)/len(query))` and its
  inverse. Under a global NW alignment that is a reasonable proxy, but the names
  promise alignment-derived coverage and do not deliver it. It is not what failed
  here — identity was — but it should be renamed or computed properly.

### Next steps

1. **Confirm the `CAGCTT` question** above, then apply the motif change (or the
   reference change) and re-baseline. This is the single highest-value open item:
   it more than doubles QC yield.
2. **Wire `coding_start`/`coding_end` into configuration** so reading-frame and
   internal-stop QC actually run. Until then, no result is checked for frameshifts
   or premature stops.
3. Re-check `low_depth` after the full run rather than on a 0.92% sample.
4. Rename or fix the coverage metrics in `qc.py`.

---

## 2026-08-12 — Standalone repository, k-mer prefilter, and a reproducibility fix

### Context inherited at the start of this session

The repository existed as `Nanopore2/Nanopore3/`, nested inside the legacy
analysis folder. The previous session had recovered the tuned demultiplexing
parameters for the 20260506 RP experiment and rebuilt the demultiplexing stage,
reducing it from 299 s to 55 s on a 20,000-read pilot. It flagged reference
assignment as the next bottleneck at 294 s and recommended optimising it before
running the full 2.17 M-read dataset.

### What changed

**1. The repository is now standalone (administrative).**
Moved `Nanopore2/Nanopore3` to `/Users/thomasfryer/Coding/Nanopore3`. It is a
self-contained Git repository; Nanopore2 and all source data were not modified.
A `.venv` was created and the package installed with `pip install -e ".[dev,report]"`.

**2. `canonical_unique_kmers` rewritten (performance).**
The old implementation called `normalize_sequence` and `reverse_complement` once
per sliding window, so a 1,291 nt read cost 4.24 ms. The new one normalises
once, splits on ambiguity codes, and reverse-complements each unambiguous run a
single time, indexing into it per window: **4.24 ms → 0.47 ms, a 9× speedup**.
Verified bit-identical to the previous implementation across 400 randomised
sequences covering ambiguity codes, empty input, and k from 1 to 15.

**3. Reference assignment restructured around the k-mer prefilter (mixed).**
The old path called `assign_sequence` once per k-mer size, and each call set
`fallback_align_all=True`. When a shortlist failed the identity/coverage floors
it aligned the read against **every** reference in the library, and then the
outer loop repeated that whole-library sweep for k=11 and k=9. Profiling
measured **78 alignments per read** on the pilot, of which 21,150 out of 23,550
came from the brute-force sweep.

The new `assign_read()` in `src/nanopore3/assignment.py`:

- extracts the insert **once** per read instead of once per k-mer size;
- caches alignment geometry per read, so each candidate reference is aligned at
  most once no matter how many k-mer sizes are consulted;
- skips a k-mer size that proposes nothing the larger sizes did not, because
  alignment evidence does not depend on k;
- replaces the whole-library sweep with a bounded **k-mer rescue**: the smallest
  k re-shortlists with no score floor, taking up to `rescue_candidates`
  (default 25) references that share any k-mer at all.

The result is **4.1 alignments per read instead of 78**. `rescue_policy` is
configurable per reference library — `kmer` (default), `all` (the old
exhaustive sweep), or `none`.

**4. `ReferenceIndex.shortlist` made hash-seed independent (scientific
correctness — this was a real bug).**
`shortlist` accumulated `scores[aliases] += 1.0 / len(owners)` while iterating a
`frozenset` of k-mers. Set iteration order for strings depends on
`PYTHONHASHSEED`, which is randomised **per process**. Floating-point addition
is not associative, so the same read scored differently in different worker
processes. Scores now accumulate as exact integer tallies per specificity class
and are summed in a canonical order.

This was not hypothetical: the first process-backend pilot differed from the
thread-backend pilot on 2 of 15,126 reads (`second_identity` and
`identity_margin`). No final call changed in that instance, but the mechanism can
move a candidate across the `minimum_kmer_score` cutoff or reorder the
shortlist, so it could change calls on other data.

**5. Parallel backend corrected to processes (performance).**
See "Lessons" — threads were actively harmful. `configs/runs/*.yaml` now use
`backend: process, jobs: 8`. `_ordered_map` gained an `initializer`/`initargs`
path so process workers build their own reference indexes once instead of
receiving them with every task, and the assignment stage now submits coarse
32-read batches rather than one future per read.

**6. Zero-length basecalls no longer abort the whole file (scientific, opt-in).**
Validating the full production profile failed outright:

```
error: .../None_sample_1.fastq: FASTQ record 91708: sequence is empty
```

Inspecting the file shows a structurally valid four-line record whose sequence
and quality are both zero length — a normal basecaller artefact, not corruption.
The strict reader treated it as fatal, which blocked the entire 2.17 M-read run.

A zero-length read is nonetheless indistinguishable from a truncated file at the
point of parsing, so the default remains a hard error and **nothing changes
silently**. `library.allow_empty_reads: true` opts in, and such reads are then
kept and classified by demultiplexing as their own counted `empty_read` state
rather than being dropped or mistaken for a barcode failure. It is enabled in
`configs/runs/20260506_lab_biotin_r1.yaml` only. With it set, the production
profile validates: **2,167,558 records**.

A direct scan of the source FASTQ counts **91 zero-length reads in 2,167,558**
(0.004%), the first at record 91,708:

```bash
awk 'NR%4==2 && length($0)==0 {n++} END {print n+0}' None_sample_1.fastq
```

The yield impact is negligible; the availability impact was total, since one such
record aborted the entire run. Confirm the `empty_read` count in the
demultiplexing summary matches 91 after the full run — a larger number would mean
the reader is now accepting something it should not.

**7. New tests and tooling.**
`tests/test_assignment_determinism.py` (5 tests) covers hash-seed stability,
shortlist order-independence, k-mer rescue versus exhaustive alignment, and the
cascade's reported k. `tests/test_empty_reads.py` (4 tests) covers the empty-read
policy. `scripts/benchmark_assignment.py` mirrors the existing demux benchmark.
`configs/runs/20260506_lab_biotin_pilot20k.yaml` makes the 20k pilot reproducible
from a committed file. Suite is 38 tests, all passing.

### Measurements

Machine: MacBook Pro, Apple Silicon, 16 logical CPUs, macOS 26.6.1, Python
3.12.10, edlib 1.3.9.post1. Workload: `runs/pilot_inputs/20260506_lab_biotin_byte_sample_20000.fastq`
(20,000 reads; 15,126 survive demultiplexing).

Per-stage wall time, from the `created_utc`/`completed_utc` fields of each
`stages/*/manifest.json`:

| Stage | Before (legacy-optimised) | After (process, 8 workers) |
| --- | --- | --- |
| 02_demux | 55 s | 23 s |
| 03_assignment | 294 s | 15 s |
| Total | 353 s | 42–54 s |

**Reference assignment is ~19× faster.** Total pipeline wall time varies between
42 s and 54 s across repeats on a loaded laptop; the stage benchmarks below are
the more reliable numbers.

Backend sweeps, identical calls in every case:

```
scripts/benchmark_demux.py       --chunk-reads 1000
  serial          591 reads/s
  thread, 4       396 reads/s     <- slower than serial
  process, 4    1,517 reads/s
  process, 8    1,818 reads/s
  process, 12   1,784 reads/s

scripts/benchmark_assignment.py  (steady state, startup cost removed by
                                  solving T(n) = startup + n/rate at n=3,000
                                  and n=15,126)
  serial          660 reads/s
  process, 4    1,378 reads/s
  process, 8    2,127 reads/s
  process, 16   2,021 reads/s     <- no gain, ~16 s of startup
```

Process workers cost roughly 1 s each to start (spawn, import, build indexes).

### Equivalence evidence

Compared against `runs/20260506-lab-biotin-pilot-20k-legacy-optimised`:

- `demux_calls.csv.gz`, `consensus.fasta`, `consensus.csv.gz`, and `qc.csv.gz`
  are **byte-identical**.
- All 15,126 assignment statuses and `reference_ids` are identical. Stage
  summaries match exactly: 12,145 `assigned_unique`, 1,901 `no_match`, 973
  `motif_missing`, 107 `ambiguous`.
- `rescue: kmer` and `rescue: all` were run over all 15,126 reads and agreed on
  every read, at 722 vs 220 reads/s serial.
- `assignment_calls.csv.gz` does differ from legacy, in recorded evidence only.
  All 740 rows with a changed `best_identity` are `no_match`. On assigned rows,
  `best_identity` changed in neither direction (0 higher, 0 lower) and
  `identity_margin` is completely unchanged — so the k-mer shortlist never
  missed a reference the exhaustive sweep would have chosen, and the ambiguity
  decision rests on the same evidence. 3 `no_match` rows carry a more accurate
  reason code (`no reference met the k-mer evidence threshold` instead of
  `best alignment did not meet identity and coverage thresholds`).
- The serial and 8-process runs are byte-identical on **every** artifact,
  including `assignment_calls.csv.gz`
  (`runs/20260506-pilot20k-fixed-serial` vs
  `runs/20260506-pilot20k-fixed-process8`).

### Lessons

**Threads cannot speed up this pipeline; processes can.** A synthetic
microbenchmark shows `edlib.align` releasing the GIL and scaling 3.8× on 4
threads, which is presumably what motivated the earlier `backend: thread`
default. That benchmark does not represent the real workload. Once the
brute-force alignment sweep was removed, edlib became a *minority* of assignment
runtime — the remainder is Python-level k-mer and dictionary work that holds the
GIL. Measured on the real workload, 4 threads made demultiplexing **33% slower
than serial** (396 vs 591 reads/s). The earlier note that "the measured thread
knee is around 2–4 workers" was measuring thread overhead, not a hardware limit.
Benchmark the real stage, never a proxy.

**`additionalEqualities` is not a cost.** Passing the 160-pair IUPAC equality
table to every edlib call was a suspected GIL bottleneck. It is not: with and
without it, runtime and thread scaling are indistinguishable. Ruled out by
measurement rather than by assumption.

**A brute-force fallback hides behind a good prefilter.** The k-mer prefilter was
already implemented and correct. Nearly all the cost came from what happened when
it *declined* to shortlist anything, compounded by retrying that fallback once
per k-mer size. When optimising, look at the failure path, not just the happy path.

**Set iteration order is a reproducibility hazard whenever floats are
accumulated.** Anything of the form `total += weight` inside `for x in some_set`
is process-dependent. This class of bug is invisible in single-process runs and
appears only once multiprocessing is enabled — precisely when it is hardest to
attribute.

**Strictness has to distinguish "malformed" from "unwelcome".** Rejecting a
zero-length read as a format error stopped a 2.17 M-read run over one record that
was, structurally, perfectly valid FASTQ. Strict parsing is right, but the
strictness belongs where the policy is visible and configurable, not buried in
the reader. The compromise is to keep the loud default and make the alternative
an explicit, counted state.

**A regression test must be shown to fail.** The first version of the hash-seed
test passed against the *unfixed* code because it rounded scores to 12 decimal
places, hiding last-place differences. It was only trustworthy after being run
against the reverted implementation and observed to fail. Always verify a
regression test detects the bug it claims to cover.

### Next steps

1. **Run at full scale.** `configs/runs/20260506_lab_biotin_r1.yaml` now
   validates against the full 2,167,558-read FASTQ and has not yet been run.
   Projected from the measured rates: roughly 20 min demultiplexing and 13 min
   assignment, versus an estimated 10+ hours before this session. Note that
   preflight alone takes about 5 min, because it hashes and scans the whole file.
   Confirm memory stays bounded with 8 workers — each holds its own reference
   indexes (~3 libraries, 1,175 references). The demultiplexing summary should
   report exactly 91 `empty_read` calls.
2. **Investigate the 973 `motif_missing` and 1,901 `no_match` reads.** That is
   19% of demultiplexed reads discarded. Determine whether they are genuine
   off-target/chimeric molecules or a motif/threshold artefact before treating
   any yield number as final.
3. ~~**Consensus is the next optimisation target**: 1,731 of 2,483 groups fail as
   `low_depth` and only 286 pass QC. That is a scientific yield question (reads
   per well) more than a speed one.~~ **Wrong on both counts — corrected by the
   2026-08-12 (later) entry.** `low_depth` is an artefact of sampling 0.92% of the
   dataset, and the QC failures are a constant 6-base boundary offset, not depth.
4. **`_ordered_map` still starves workers by design.** It keeps exactly `jobs`
   futures in flight. Increasing the queue depth was measured and made no
   difference at current batch sizes, but it will matter if per-batch cost drops.
5. **Address the pre-existing lint debt.** `ruff check src/` reports 52 findings,
   almost all pre-existing; CI does not currently run ruff. Either gate CI on it
   or record the decision not to.
6. **Consider a `nanopore3 benchmark` subcommand** so the two benchmark scripts
   are discoverable, and add a golden serial-vs-process test to CI to protect the
   determinism fix at the pipeline level, not just the unit level.

### Open questions for the experiment owner

- `configs/runs/20260506_lab_biotin_r1.yaml` maps plates RP01–RP07 to reference
  libraries, inferred from legacy output filenames with only RP07 confirmed. The
  remaining six mappings still need confirmation before results are reported.
- The pilot input is a 20,000-read byte-offset sample, not a random sample. It is
  fine for performance work and for equivalence checking, but it is not a
  statistically representative subset.
