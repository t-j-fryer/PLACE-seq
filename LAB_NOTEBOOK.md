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

## 2026-08-21 (third) — v7: insert-scoped acceptance floors, and the artefact is gone

`runs/260608-full-length-v7`, 18.6 min - the same as v6b, so the extra alignment
per candidate costs nothing measurable.

### What changed

**Identity and coverage floors now apply to the insert region, not the whole
amplicon.** A full-length reference is ~70% constant flank, so a whole-amplicon
metric is roughly 3.5x less sensitive to anything wrong with the insert than the
same floor was when the reference *was* the insert. A read carrying 250 nt of a
different design scored 0.927 identity and 0.862 coverage over the amplicon,
clearing floors of 0.80 and 0.70, and then built a designed consensus that
reported a perfect match.

The library's own floors are re-applied where they discriminate. A read that
fails gets `assignment_status: insert_mismatch`, is kept out of that design's
consensus, and carries `insert_identity` and `insert_query_coverage` in
`assignment_calls.csv.gz` so the decision is auditable. `insert_thresholds:
false` opts a library out.

Chimeric-read exclusion is also on (`exclude_reads: written`).

### It removes the artefact and leaves everything else alone

| | v5c (neither) | v6b (exclusion) | v7 (both) |
|---|---|---|---|
| reads held back at assignment | - | - | **15,635** (2.5%) |
| chimeric clones written | 90 | 103 | **103, byte-identical to v6b** |
| designed consensuses | 3,463 | 3,406 | **3,334** |
| `mixed_variants` | 128 | 108 | **52** |
| QC failures | - | 366 | **301** |
| consensuses whose reads leave >2% unexplained | 7 | 1 | **1** |

Against v6b: **72 removed, 0 newly appearing, 98.2% of survivors byte-identical.**
The 72 had median identity 86.4% to their supposed design and 69 of 72 were below
95%.

**The chimera stage summary is byte-identical between v6b and v7**, which is the
check that mattered: the gate operates on assignment, and chimera detection reads
the demultiplexed reads directly, so it must be unaffected. It is.

`reads_excluded_as_chimeric` fell from 7,849 to 4,806, because the gate now
catches most of those reads earlier. Exclusion has become belt-and-braces, as
predicted.

### The three high-identity removals, checked one at a time

Having got this wrong once by not checking, all three were examined:

| well | design | v6b | insert-scoped | verdict |
|---|---|---|---|---|
| RP05 D9 | `B_Block_1_dTF085_84_3` | 99.1% from 7 reads | median coverage **0.66** | correctly removed - the reads do not cover this design |
| RP05 F4 | same design | 99.9% from 7 reads | coverage 1.00, identity 93.7%, 2/7 reads held | **marginal clone lost to the depth floor** |
| RP05 E7 | same design | 99.4% from 7 reads | coverage 1.00, identity 93.8%, 2/7 reads held | same |

So the cost is **two marginal clones out of 3,406**, each supported by 7 reads,
which fell to 5 when two of their reads were individually held back - below
`minimum_depth: 6`. That is a depth effect, not a misjudgement: the reads that
went were the ones that individually failed the floor. Recoverable by lowering
`minimum_depth` if those clones matter.

### Independent validation: the replicate concordance holds

Two culture plates sequenced twice, once on their own barcode and once inside the
pool:

| | RP06 / SUMO_A_P1 | RP07 / SUMO_A_P2 |
|---|---|---|
| dedicated wells matched | 82/82 | 93/93 |
| clones byte-identical over ~1.2 kb | 83/84 | **96/96** |

**179 of 180 clones byte-identical from two independent read sets.** The one
exception is RP06 F08, 9 edits - a well that also holds a chimeric clone, so the
gate held back slightly different reads on each barcode. This is the check that
the gate did not quietly corrupt anything, and it passes.

Design recovery is unchanged - 75.0 / 86.0 / 80.4 / 96.5% perfect, the same as
v5c and v6b - and so is per-region accuracy: 99.997% / 99.991% / 99.999% with
96.4% of QC-passing clones base-perfect across the whole amplicon.

### A result that came free

Apparent polyclonality drops, because some of it was the artefact:

| Library 3 (B), distinct sequences per well | v5c | v7 |
|---|---|---|
| exactly one | 84.3% | **88.2%** |
| two | 11.8% | **7.7%** |
| three or more | 1.1% | **0.6%** |

Wells that looked like they held two clones were in part holding one clone and
one description of a chimera. The polyclonality estimate is now cleaner.

### Figures

All four regenerated into `runs/260608-full-length-v7/figures/`.

### Lesson

**A threshold is only meaningful relative to the sequence it is measured over.**
Moving from 350 nt inserts to 1,200 nt amplicons kept every number in the config
and silently weakened all of them by 3.5x. Nothing failed, nothing warned, and
the symptom appeared three stages later as a clone that did not exist. When the
unit of measurement changes, every threshold expressed in it needs re-deriving -
not just the ones that look related.

### Next steps

1. **Runs: keep `v7`.** `v5c` and `v6b` are superseded - say the word and they go.
2. The 132 PCR-origin chimeras are counted but not localised.
3. Scaffold synthesis fits its own clone within 10 edits in only 10 of 103 cases;
   affects the chimera consensus, not whether a molecule is called chimeric.
4. RP01-RP04 and RP09-RP11 still have no full-length config.

---

## 2026-08-21 (second) — CORRECTION: the chimera calls were right, my analysis was not

This corrects the 2026-08-21 entry, which concluded that read exclusion deleted
four genuine clones and should ship disabled. **Both conclusions were wrong.**
`exclude_reads` is back on.

### What I got wrong, and how

I read "146 reads, 100% identity" as "146 reads that perfectly match this
design". It is not. It is *one consensus*, built from 30 of those reads, that
aligns to the design with zero edits. Those are different claims, and the
difference is the whole problem.

Checking each of the four one at a time:

| well | design length | read length (median) | dominant signature | parents' identity to each other |
|---|---|---|---|---|
| RP05 A10 | 285 nt | 458 nt | 145/146 reads, 2 segments | **52%** |
| RP05 E6 | 321 nt | 588 nt | 112/113 reads, 2 segments | **53%** |
| RP05 H1 | 405 nt | 508 nt | 244/246 reads, 2 segments | **48%** |
| RP08 C6 | 330 nt | 526 nt | 76/79 reads, 2 segments | **48%** |

Every one is a chimera of two *unrelated* designs, in reads 150-250 nt longer
than the design they were assigned to, agreed by essentially every read in the
well. The chimera detector was right in all four cases. **The exclusion did not
delete genuine clones; it deleted four chimeras that were being reported as
perfect designed clones.**

I also tested the detector against pristine input, which I should have done first:
profiling all **1,016 error-free reference sequences** through it, every single one
reads as exactly itself. Zero false positives. My "similar designs confuse the
window calls" story had no basis - and it was checkable in one command.

The "50 of 103 calls explained as well by a single design" measurement stands as
a number but not as an interpretation: it measures **scaffold quality**, not
whether a molecule is chimeric. The acceptance test I proposed on the back of it
- reject a chimera unless its spliced pair beats the best single reference -
would have **rejected these four true chimeras**, causing exactly the error I
thought I was preventing. Dropped.

### Why a chimera reports 100% identity to one of its parents

This is the finding worth keeping, and it is not about chimeras:

**A reference-guided consensus represents the part of a molecule that aligns to
the reference, and silently drops the rest.** The consensus takes its length from
the reference, so 250 nt of foreign sequence in the reads has nowhere to appear.
Identity comes out 1.000. Coverage cannot see it either: QC compares consensus
length to reference length, and those are equal by construction.

And the thresholds that should have stopped the reads did not, because of
full-length references:

| metric, over the whole amplicon | median for those 146 chimeric reads | floor |
|---|---|---|
| identity | 0.927 | 0.80 |
| query coverage | 0.862 | 0.70 |
| reference coverage | 1.000 | 0.70 |

The constant flanks are 849 of ~1,200 bases - **70% of every reference is
identical across the whole library**. A 250 nt foreign segment is 17% of the
molecule and cannot pull a whole-amplicon metric below a floor that was set when
the reference *was* the insert. Over the insert alone that read covers ~0.53 and
would fail outright.

### The clean fix, in the pipeline's own terms

**Apply the identity and coverage floors to the insert region, not the whole
amplicon.** The pipeline already locates the insert span - that is how the region
metrics work. Then a read carrying foreign insert sequence fails at assignment on
its own merits: no cross-stage plumbing, no read claiming, and it cannot cost a
genuine clone, because a genuine clone's insert covers its design.

That makes read exclusion belt-and-braces rather than the only defence. Not
implemented: it changes what gets assigned, so it needs its own run and its own
comparison. **It is the right next change.**

### Measured effect of exclusion, now that it is trusted

Reads whose alignment leaves part of the molecule unexplained, per designed
consensus:

| | v5c (no exclusion) | v6b (exclusion on) |
|---|---|---|
| consensuses with >2% unexplained | 7 | **1** |
| consensuses with >10% unexplained | 5 | **0** |
| worst | 0.773 coverage | 0.932 |

### Lessons

1. **"100% identity" is a statement about a consensus, not about a molecule.** A
   reference-guided consensus cannot report what it has no coordinates for.
2. **Test a detector on perfect input before theorising about why it fails.** One
   command over 1,016 references would have killed my false-positive story
   immediately.
3. **When thresholds move to a longer reference, they get weaker.** Moving from
   350 nt inserts to 1,200 nt amplicons diluted every insert-level defect by 3.5x
   against floors nobody re-derived. Which is the actual bug.

### Next steps

1. **Insert-scoped identity and coverage floors at assignment** - the fix above.
2. Then re-check whether `exclude_reads` is still needed at all, or has become
   redundant.
3. Scaffold synthesis fits its own clone within 10 edits in only 10 of 103 cases;
   worth improving, but it affects the chimera consensus, not whether a molecule
   is called chimeric.

---

## 2026-08-21 — The chimeric-read exclusion works, and must not be switched on yet

Implemented as agreed, measured, and then **defaulted off**, because measuring it
exposed that the chimera calls it rests on are not sound enough to rest on.

### What was built

1. **Chimera detection moved ahead of consensus** (`04b_chimera` -> `03b_chimera`).
   It only ever needed the demultiplexed reads and the reference list, so nothing
   is lost by running it first.
2. **It records the reads it claims** (`clones/claimed_reads.csv`: read_uid,
   chimera_id).
3. **`chimera.exclude_reads`** = `written` | `all` | `none`. With `written`, a read
   belonging to a written chimeric clone stops also voting for a parent design -
   it is not discarded, it is already building its own clone's consensus.
4. **The consensus stage reports what it skipped** (`reads_excluded_as_chimeric`)
   rather than skipping silently, and its fingerprint includes both the setting
   and the digest of the claimed-read list.
5. **Profiling recruits more reads.** The insert boundaries are searched first and
   only if that fails is the insert sliced positionally out of the primer-anchored
   region the assignment already matched. That recovers reads that were assignable
   but not profilable - the six in RP07 A12 that escaped detection and built a
   spurious consensus.

### The exclusion does what it was meant to

`runs/260608-full-length-v6b` against `runs/260608-full-length-v5c`:

| | v5c | v6b |
|---|---|---|
| chimeric clones written | 90 | **103** (0 lost, 13 new) |
| designed consensuses | 3,463 | 3,406 (**57 removed, 0 appearing**) |
| of the shared, byte-identical | - | **3,394 / 3,406 (99.65%)** |
| reads excluded as chimeric | - | 7,849 |

The 57 removed had median identity 92.9% to their supposed design, 47 of them
below 95%, 26 already graded `mixed_variants`. Exactly the population it was
aimed at.

### But it also deleted four genuine clones

Four of the 57 were at 98-100% identity, and not marginal: RP05 A10 held **146
reads giving a 100%-identity consensus** of `B_Block_2_dTF085_DENOISE_263_4`, and
it disappeared, because a 160-read "chimeric clone" in the same well listed that
design as its second parent and claimed the reads.

146 reads agreeing perfectly with one design **are** that design. So that chimera
call is a false positive.

### Which led to measuring the calls themselves

For each of the 103 written chimeric clones, how well does its spliced pair of
references explain its consensus, against the best *single* reference - both by
global alignment, so neither model is excused for unmatched ends:

| | |
|---|---|
| median edits saved by two references instead of one | **4** |
| clones a single design explains as well or better | **50 of 103** |
| clones fitting their own scaffold within 10 edits | **10 of 103** |
| median clone-to-scaffold distance | **87 edits** |

The good calls are unmistakable - scaffold fit of 1-20 edits, beating the best
single reference by 71-97. The bad ones are equally unmistakable: scaffold fit of
200-300 edits while a single design sits 90-100 away. Roughly half the calls are
in the second group.

**This is pre-existing** - the same 90 clones are in v5c - but it was harmless
there: a spurious chimera FASTA sat beside the true clone. Exclusion makes it
consequential, because a false chimera now deletes a real clone.

So `exclude_reads` ships as **`none`**, and the config says why. The plumbing is
in place and one word turns it on.

**I nearly reported the wrong number here.** My first measurement compared the
scaffold by global alignment against the single reference by infix alignment,
which lets the single reference ignore unmatched ends. It said 57 of 103. The
fair comparison says 50 - the conclusion held, but the number would have been
wrong and the method indefensible. *Comparing two models means aligning them the
same way.*

### The fix that unlocks this

A chimera should have to earn its call: its spliced scaffold must explain the
clone better than the best single reference, by a margin. That is one extra
alignment per candidate group in a stage that already aligns every window, so the
cost is negligible. It would drop ~50 of 103 calls, and make exclusion safe -
the four genuine clones were deleted by exactly the calls this test rejects.

Not implemented, because it changes what the pipeline calls a chimera, which is a
scientific definition and wants sign-off rather than my judgement.

### Next steps

1. **Decide the acceptance test above**, then turn `exclude_reads` on.
2. Runs kept: `260608-full-length-v5c` (the published figures) and
   `260608-full-length-v6b` (evidence for this entry). `v6` was deleted - it
   carried a worse extractor that lost 11 chimeric clones, found by the same
   comparison and fixed before v6b.
3. Carried over: RP01-RP04 and RP09-RP11 have no full-length config; PCR-origin
   chimeras are counted but not localised.

---

## 2026-08-20 (eleventh) — What is actually in RP07 A12, and one real defect

Runs pruned to `260608-full-length-v5c` alone: v4, v5 and v5b deleted, 6.5 GB
recovered. `20260506-*` and `260608-AI-DBTL-v2/v3` remain and are a separate
question - see next steps.

### The question

The one non-identical clone in the full-length replicate concordance was RP07
A12, graded `mixed_variants` in both runs. Was the k-mer assignment wrong, and
does the well hold two clones that would each make a clean consensus?

### Answer: one molecule, and it is chimeric

Added `scripts/dissect_well.py`, which clusters a well's reads **without
reference to any design** and then asks what each group matches - deliberately
the reverse of the question that produced the confusing answer.

```
python scripts/dissect_well.py --run-dir runs/260608-full-length-v5c \
    --plate RP07 --well A12
```

| evidence | value |
|---|---|
| reads in the well | 367 (293 with a usable insert region) |
| assignment status | **294 ambiguous**, 29 assigned_unique, 43 motif_missing, 1 no_match |
| reference-free clusters | **1**, zero co-varying positions, mean mismatch to centre 0.0 |
| reads sharing one positional signature | **285 of 293** |
| that signature | `A_Block_2_dTF079_0_..._50_1` >> `A_Block_2_dTF079_2_..._42_4` |
| its consensus vs the closest single reference | **77.5%** (56 edits over 249 nt) |

So:

1. **The assignment was not wrong.** It called 294 of 367 reads `ambiguous`,
   which is the correct answer when no single reference fits - and none does.
2. **There is one molecule in the well, not two.** The reads agree with each
   other: zero co-varying positions and a mean mismatch to the cluster centre of
   0.0. A mixed culture would show co-varying positions.
3. **That molecule is a chimera of two Block_2 designs.** Same block means the
   same colony, so it is a real clone, not a PCR artefact - and the chimera stage
   already found it, with a consensus byte-identical between RP07 and RP05.
4. **Splitting the reads would not give two clean clones.** Each parent covers
   only part of the molecule; the correct object is the spliced one, which the
   pipeline already writes.

### The real defect this exposed

The 19 reads assigned to `A_Block_2_dTF079_2_..._42_4` at **89.7% identity**
passed the library's `minimum_identity: 0.80` floor and went on to build a
designed-clone consensus - the `mixed_variants` one, 67 insert edits from its
supposed design. **That consensus should not exist.** It is the residue of a
chimeric clone described as a design.

Scale, across v5c's 3,463 designed consensuses:

| | count | share |
|---|---|---|
| identity <95% to their design | 189 | 5.5% |
| of those, in a well that also produced a chimera | 118 | 62% |
| baseline: clean consensuses in a chimera-bearing well | 1,333 | 40.7% |

A 1.5x enrichment - real but modest, and the pooled RP05 wells make the baseline
high because almost every well holds many clones. So this is a ~5% tail, not a
systemic problem, and `mixed_variants` is already flagging most of it (79 of 189).

**The fix worth making:** the chimera stage knows which reads belong to a written
chimera group, and the consensus stage does not. Passing that set across would let
a designed-clone consensus either exclude those reads or be flagged as
chimera-derived. Not done - it changes consensus membership and so needs its own
run and comparison.

### Two mistakes I made getting here, both caught by checking

1. **I wrote my own window scan and had the alignment backwards** - aligning each
   whole 249-351 nt reference *into* a 90 nt window, which is meaningless. It
   reported a five-segment chimera of Block_1 designs, contradicting the chimera
   stage's two Block_2 parents. The contradiction is what exposed it. Replaced
   with `chimera.signature_for`, the tested profiler the chimera stage itself
   runs, which gave the two-segment answer that agrees.
2. **I assumed the differing clone in that well was the chimera.** It is the
   opposite: the chimera is byte-identical between the two runs and the *designed*
   clone is the one that differs. Checked before writing it down.

**Lesson: when a diagnostic disagrees with a pipeline stage that answers the same
question, suspect the diagnostic.** The stage has tests.

### Next steps

1. **Decide on the chimera-read exclusion above.** It would remove ~118 spurious
   designed consensuses and cost one run.
2. **`260608-AI-DBTL-v2` and `v3` (4.2 GB) and the `20260506-*` pilot runs (2.2
   GB) are still on disk.** v2 and v3 are superseded runs of this dataset; the
   20260506 set is different data. Say which to drop.
3. Carried over: RP01-RP04 and RP09-RP11 have no full-length config; the 72
   PCR-origin chimeras are counted but not localised.

---

## 2026-08-20 (tenth) — v5c: RP06 and RP07 join the full-length run

`runs/260608-full-length-v5c`, 21.3 min. `configs/runs/260608_full_length.yaml`
(renamed from `260608_rp05_rp08_full_length.yaml`, which the barcode list had
outgrown) now maps RP05, RP06, RP07 and RP08 - all four share this construct.
RP06 and RP07 are dedicated barcodes for two of RP05's pooled culture plates, so
including them makes the concordance check run over the whole amplicon.

### The strongest validation in the set

Each of those two culture plates was sequenced twice, independently: once on its
own barcode and once inside the 22-plate RP05 pool.

| | RP06 / SUMO_A_P1 | RP07 / SUMO_A_P2 |
|---|---|---|
| dedicated wells matched | 82/82 (100%) | 93/93 (100%) |
| shared clones byte-identical over ~1.2 kb | **84/84 (100%)** | **96/97 (99.0%)** |
| wells recovered only from the pool | 13 | 2 |
| wells recovered only from the dedicated barcode | 0 | 0 |

**180 of 181 clones are byte-identical across the full amplicon from two
independent sets of reads** - roughly 216 kb of independently reconstructed
sequence with one disagreement. At the insert level (v4) plate 1 had one
non-identical clone; over 1.2 kb it has none.

### The one disagreement is a well the pipeline already flags

RP07 A12 / `A_Block_2_dTF079_2_SUMO_l117_s583407_mpnn1_model2_42_4`, 32 nt apart
between the two runs:

- both runs grade it **`mixed_variants`** - the reads for this design disagree
  with each other beyond the support threshold, so the pipeline says in both runs
  that this is not a single clean clone;
- both carry ambiguity codes (24 bases from 12 reads in the pool, 8 from 19 on
  the dedicated barcode) and both differ from the design (78 and 67 insert
  edits);
- **all 32 differences are inside the insert. Both constant regions are
  byte-identical in both runs, and both match the reference exactly.**

That last line is what full-length references buy. Two independent samples of a
heterogeneous well give two different consensuses over the variable region while
agreeing perfectly over the 849 constant bases - so the disagreement is the well,
not the method, and the reader can see that rather than take it on trust. The
same well also holds a chimeric clone, which *is* byte-identical between the two
runs.

### Figures

All four, in `runs/260608-full-length-v5c/figures/`:
`replicate_concordance`, `platform_sequences_per_well`,
`platform_reference_recovery`, `platform_sequence_populations`, plus their CSV and
JSON companions and `region_summary.json`. The report stage's `fig1`-`fig4` are in
`stages/06_report/figures/` as always.

Numbers unchanged from v5b where they should be: insert-scoped recovery 75.0 /
86.0 / 80.4 / 96.5 and per-region accuracy 99.997% / 99.991% / 99.999% over 3,105
QC-passing consensuses, 96.4% base-perfect across the whole amplicon. Adding two
barcodes adds clones without moving the per-clone measurements, which is the
expected result and worth having checked.

### Next steps

1. **RP01-RP04 and RP09-RP11 still have no full-length config.** RP01-RP04 are
   opTF001 on the same 5' motif but a different insert set; RP09-RP11 end at a
   different 3' constant region and would each need their own flanks.
2. The 72 PCR-origin chimeras are counted but not localised (carried over).
3. **v4, v5, v5b and v5c now coexist**, 2-3 GB each. v5c supersedes v5b, which
   superseded v5. Worth pruning once the figures are settled - say which to keep.

---

## 2026-08-20 (ninth) — v5b figures, and a like-for-like fix

Figures for `runs/260608-RP05-RP08-v5b`. Two sets, in two places:

**From the run itself**, `stages/06_report/figures/`: `fig1_plate_occupancy`,
`fig2_clonality`, `fig3_deconvolution`, `fig4_culture_plates`. These are written
by the report stage and were there as soon as the run finished.

**From the analysis scripts**, `figures/`: `platform_sequences_per_well`,
`platform_reference_recovery`, `platform_sequence_populations` (PDF, 600-dpi PNG,
transparent SVG each), plus `platform_reference_recovery.csv`,
`platform_comparison.json` and `region_summary.json`. These had not been
regenerated - I had raised it as an open question rather than doing it.

`replicate_concordance` **cannot** be made for v5b: it needs RP06 and RP07, and
this config maps only RP05 and RP08. Extending it means adding those barcodes to
the full-length config, which is a config change and a fresh run.

### The fix generating them forced

Scoring the nanopore side over the whole amplicon while Illumina only ever reads
the insert is not a comparison. The v5b figures first came out with nanopore held
to a standard three times longer than its opponent:

| | scored on the amplicon | scored on the insert |
|---|---|---|
| Library 1 (stuffer) Perfect recovery | 73.5% | **75.0%** |
| Library 3 (A) | 85.1% | **86.0%** |
| Library 3 (B) | 80.4% | **80.4%** |
| Library 3 (A+B) | 96.2% | **96.5%** |
| stuffer consensus Perfect population | 88.8% | **91.5%** |
| Library 3 (A) population | 87.4% | **90.2%** |

`compare_platforms.py` now takes `--scope insert|amplicon`, defaulting to
`insert`, and a full-length run supplies the insert-only outcome per consensus
from the region columns: perfect when the insert matches exactly, screenable when
it does not but frame and stop checks pass, otherwise other. An insert-only run
is unaffected - it has nothing else to score.

Insert-scoped v5b now reproduces v4 to within 0.3 points on recovery (75.0 /
86.0 / 80.4 / 96.5 against 75.0 / 85.7 / 80.4 / 96.5), which is the right answer:
the two runs demultiplex and deconvolve identically, so a figure about *which
designs were recovered* should not move when the consensus gets longer.

The population figure is 2-3 points lower than v4 for a real reason, not a
scoping one: v5b includes the 188 wrecked clones v4 never assembled, and they are
observations, so they belong in a population.

**Lesson: when the measurement changes length, check every figure that compares
it to something else.** The recovery numbers moving by ~1.5 points was the only
visible symptom, and it would have read as noise.

### Also worth knowing

`platform_sequences_per_well` moved for a real reason too - v5b finds a clone in
more wells (empty wells 1.6% against 3.5% for the stuffer library) and finds more
polyclonal wells (Library 3 (B) two-sequence wells 11.8% against 6.1%), both
because the primer anchors are present in 92% of reads where the insert motifs
were in ~85%. More reads per well crosses `minimum_depth` more often, for second
clones as well as first.

### Next steps

1. **Decide whether RP06/RP07 join the full-length config**, which would restore
   the replicate concordance figure at full length - the strongest validation in
   the set, since it would then compare 1.2 kb of independently sequenced
   sequence rather than 350 nt.
2. Carried over: the 72 PCR-origin chimeras are counted but not localised, and
   RP01-RP04 and RP09-RP11 have no full-length config.

---

## 2026-08-20 (eighth) — Chimera detection back on in v5, scoped to the insert

**Scientific change**, run as `runs/260608-RP05-RP08-v5b` (v5 plus chimeras;
`runs/260608-RP05-RP08-v5` is kept as the no-chimera comparison).

### Why it needed scoping rather than just switching on

Positional k-mer profiling slides a window along a read and asks which
references match it.  That question only carries information where the
references differ.  In the full-length library 849 of ~1,200 bases are identical
across all 684 references, so every window in a constant region matches
everything and contributes nothing but noise to a signature.

So detection runs on an insert-only view of the same run:

* the inserts are sliced back out of the full-length references by the same span
  the flanks define - no file is re-read, and the sequences are identical to the
  insert FASTA by construction;
* the region is bounded by the *inner* ends of the constant regions
  (`Flanks.insert_anchors()`) instead of the primer sites.

`insert_anchors()` rejects an anchor that recurs in the constant regions, since
that could bound the wrong region - the same guard the primer anchors already
had.

Clones written this way are insert-scoped while their designed neighbours span
the whole amplicon, which is a mixed convention inside one run.  Rather than
leave it to be inferred from a length, each carries `region=insert` in its FASTA
header and a `region` column in `clones.csv`.

### It agrees with the insert-only run

| | v4 (insert-only run) | v5b (full-length run, insert-scoped) |
|---|---|---|
| assembly-origin clones written, RP05+RP08 | 84 | 88 |
| PCR-origin detected, not written | - | 72 |
| same well + same parent pair | **83 of 84** | |
| **same junction window** | **83 of 83** | |
| byte-identical consensus | 64 of 83 | |

Every call the two runs share agrees on *where the junction is*, which is the
biological claim.  The 19 whose sequence differs do so because the contributing
read set differs - v5b keeps reads carrying both primer sites, v4 kept reads
carrying both insert motifs - and a different subset of a well's reads gives a
consensus a few bases different.  Both are valid consensuses of the same clone.
Median difference 5 nt; the largest, 48 nt, is a clone where the read sets
diverge most.

The 5 calls v5b makes that v4 did not all sit at 6-7 reads, right on
`minimum_depth: 6`: the full-length anchors are present in 92% of reads against
~85% for the insert motifs, so a few groups tip over the threshold.  The single
v4-only call sits at 6 reads for the same reason in reverse.

All 88 resolve to a culture plate, all are two-parent, and all are
assembly-origin - both parents in the same block, so the same colony - which is
the only kind that can be a real clone.  72 PCR-origin products were counted and
not written, as before.

### Cost

17.4 min end to end against 14.3 without chimeras, so the stage costs 3.5 min -
half of v4's 6.7, because only two barcodes are mapped here rather than eleven.

### Lesson

**A stage that depends on references being discriminative needs to be told which
part of them is.** Switching chimera detection on unchanged would have run it
over 1,200 bases of which 850 are shared, and it would not have crashed or
warned - it would have produced worse signatures and fewer calls, and the config
said nothing about why.  The scoped view is now explicit in the config comment
and in the function's docstring, and `region=` on every clone says which
convention produced it.

I also guessed wrong once here: I expected the differing consensuses to be
explained by the 6-base boundary convention between v4's reverse motif and v5b's
insert anchor, and tested it - 1 of 19. The real cause was the read sets. Worth
recording that the test took a minute and stopped a wrong sentence reaching the
notebook.

### Next steps

1. **The 72 PCR-origin products are still only counted.** In a full-length run
   their junctions could be localised across the whole amplicon rather than the
   insert; whether that is worth a second, amplicon-scoped pass is an open
   question.
2. **v5 and v5b differ only by the chimera stage.** If disk matters, `v5` can go;
   it is kept for now as the clean no-chimera comparison.
3. Carried over: whether to regenerate the insert-level figures against v5b, and
   whether to give RP06/RP07 a full-length config so replicate concordance covers
   the whole amplicon.

---

## 2026-08-20 (seventh) — v5: full-length consensus between the primer sites

**Scientific change.** The consensus now spans the whole amplicon between the two
nanopore primer binding sites and is scored against full-length references, built
by joining each insert from the oligo-pool tool to the constant regions either
side. v4 reconstructed ~350 nt of insert and discarded ~850 sequenced bases per
read; those bases are now reconstructed, and accuracy is reported per region.

Run: `runs/260608-RP05-RP08-v5`, config
`configs/runs/260608_rp05_rp08_full_length.yaml`, RP05 (SUMO A+B) and RP08 (LAB
stuffer) only - the other barcodes carry constructs whose constant regions
differ.

### How it is configured

```yaml
reference_libraries:
  sumo_ab:
    fasta: ../../runs/refs_260608/sumo_ab.fasta   # inserts, as the oPool tool emits
    flanks:
      upstream: cATAATCCGCACGCATCTGG...GCAGCTATGCAGCTT
      downstream: agtGGaTCC...ggattggcgaatgggacgc
      anchor_length: 20
```

or, for the friendlier route the user asked for, one example construct instead of
two pasted strings:

```yaml
    flanks:
      template: ../../runs/refs_260608/example_construct.fasta
```

from which the constant regions are derived by finding whichever insert the
example contains. Nothing else changes: the primer sites at the outer ends of the
constant regions become the region boundaries automatically, so no motif has to
be restated, and an explicit motif still wins if a construct needs one.

Given flanks are 135 nt (5') and 754 nt (3'), the anchors take 20 nt from each
end, so 849 constant bases are added to every reference and full-length
references run 1,038-1,299 nt against inserts of 177-450 nt.

### Checked before trusting it

- **The reads reach both primer sites.** 91.7% of RP05 and 92.8% of RP08 reads
  contain both anchors within 3 edits, median anchored span 1,227 nt against an
  expected 1,240 for a 351 nt insert.
- **The flanks are the right ones.** Raw read identity to the full-length
  references is 97.6% median - a wrong constant region would have collapsed this.
- **The built reference is what the pipeline extracts.** A round-trip test builds
  a synthetic read, runs `extract_insert` with the derived anchors and asserts the
  result equals the built reference, so the "inner flanks" convention cannot
  silently drift.
- **The ORF spine is whole codons.** 9 nt upstream of the insert plus 555 nt to
  the terminal stop = 564, divisible by 3, so frame reduces to the insert's own
  length as it did in insert mode.

### Result: QC-passing clones are now verified across the whole amplicon

2,931 of 3,284 consensuses (89%) grade `perfect` or `screenable`. For those:

| region | length | mean identity | exact | errors/kb |
|---|---|---|---|---|
| 5' constant | 115 nt | 99.998% | 99.8% | 0.02 |
| insert | 315 nt | 99.991% | 97.4% | 0.09 |
| 3' constant | 734 nt | 99.999% | 99.2% | 0.01 |
| whole amplicon | ~1.2 kb | **99.996%** | **96.4%** | - |

**96.4% of passing clones are now confirmed base-perfect across the entire 1.2 kb
amplicon, not just the insert.** That is the claim v4 could not make at all. Of
the insert-perfect clones, only 1.09% differ from the reference anywhere in the
constant region, which answers the open question from earlier today: mutations
hiding in LgBit are rare, and now they are measured rather than assumed.

The 353 clones that fail QC put their damage almost entirely in the insert:

| region | mean identity | errors/kb |
|---|---|---|
| 5' constant | 98.993% | 10.1 |
| insert | 75.110% | 257.6 |
| 3' constant | 99.180% | 8.2 |

25x the error rate in the synthesised insert against the clonal vector in the
same molecules. The vector is an internal control in every read: it says the
consensus procedure itself is near-perfect, so insert differences are real
synthesis and assembly errors rather than sequencing artefacts. That control did
not exist when only the insert was reconstructed.

### Full-length mode does not degrade the insert, and finds more clones

| | clones | insert mean identity | insert exact |
|---|---|---|---|
| called by both v4 and v5 | 3,096 | v5 **99.755%** / v4 **99.759%** | 93.6% both |
| v5 only | 188 | 57.150% | 4.8% |
| v4 only | 14 | 88.149% | 0.0% |

On every clone both runs call, the two agree to 0.004 points and 97.2% of inserts
are byte-identical (v5 has *fewer* ambiguous bases: 157 against 271). v5 finds
188 clones v4 never assembled, because the primer anchors are present in 92% of
reads where both insert motifs were present in ~85%. Those extra clones are
wrecked molecules - 57% insert identity, 4.8% exact - which is useful data about
what is in the wells, but they are not a like-for-like accuracy population.

**The mistake I made, recorded because it is the easy one to make here.** My
first summary pooled all 3,284 consensuses and reported the insert at 97.3% mean
identity and 29 errors/kb, which read as though full-length mode had made the
insert three times worse. It had not: 188 wrecked clones out of 3,284 moved the
pooled mean by 2.4 points. I only caught it by asking why two numbers over the
same clones disagreed, and confirming per clone that they did not.
`summarise_regions.py` now reports by grade so the split is visible rather than
averaged away. **A mean over a population that changed is not a comparison.**

### Performance

14.3 min end to end against v4's 20.2, on the same 2.55 M reads. Assignment cost
3.1 min against 2.8 despite references 3.5x longer, because the
specificity-weighted k-mer prefilter is unaffected: constant-region k-mers are
owned by all 684 references and weighted 1/684, so the insert k-mers still decide
the shortlist. Consensus 0.9 min against 0.6. v5 has no chimera stage, which is
where v4 spent 6.7 min.

`chimera.enabled: false` in this config on purpose: positional k-mer profiling
assumes the region discriminates between references, and 849 of ~1,200 bases are
now identical across every reference. Use the v4 insert-level run for chimera
calls.

### Added

- `src/nanopore3/flanks.py`, `configs/runs/260608_rp05_rp08_full_length.yaml`,
  `scripts/summarise_regions.py`
- `references.read_fasta` takes a sequence transform, so flanking happens before
  alias detection and digesting: the reference digest describes what is searched,
  not the file on disk. Flanking cannot merge two distinct inserts or split two
  identical ones, and there is a test for each.
- `qc.region_metrics` and `qc.coding_checks_in_consensus`; the QC table gains
  `{flank_5p,insert,flank_3p}_{length,edit_distance,identity}` columns, present
  only when a library declares flanks, so an insert-mode run keeps its schema.
- 32 tests across `tests/test_flanks.py`, `tests/test_qc_full_length.py` and
  `tests/test_full_length_wiring.py`; 241 pass.

### Next steps

1. **The insert-level figures from earlier today are still v4's.** The platform
   comparison, the replicate concordance and the culture-plate figures all read
   the v4 run. Decide whether they should be regenerated against v5 - the
   deconvolution and demultiplexing are unchanged, so only the accuracy panels
   would move.
2. **353 QC-failing clones with 258 errors/kb in the insert** are worth a look
   as a set: how many are truncations, how many chimeras, how many point
   mutations. The v4 chimera stage covers some of this; a full-length view would
   localise junctions much better than an insert-only one.
3. **RP01-RP04 and RP06, RP07, RP09-RP11 have no full-length config.** RP06 and
   RP07 share this construct and would extend the replicate concordance to the
   whole amplicon; the others need their own constant regions.

---

## 2026-08-20 (sixth) — Allocated wells set the depth; outcome populations replace mean accuracy

Three changes to the platform comparison, requested after reading the first
version. Analysis only; no stage or default touched.

### 1. Matched depth now uses allocated wells, not recovered ones

Previously Illumina was given as many reads per block as the nanopore run had
wells *yielding* a clone of that block. That charged nanopore's own failures to
Illumina as reduced depth. It now uses the wells **allocated** to each block -
the colonies actually picked.

Which well got which block is recorded nowhere, so a plate's 95 allocated wells
are apportioned among the blocks it carried in the ratio the recovered clones
show, by largest remainder, with every plate contributing exactly 95. Total
depth rises from 3,000 to **3,135 reads** (1,045 per set, over 51 blocks).

Perfect recovery, with the change:

| set | nanopore | Illumina, recovered-well depth | Illumina, allocated-well depth |
|---|---|---|---|
| Library 1 (stuffer) | 75.0% | 54.5% | **55.7%** |
| Library 3 (A) | 85.7% | 67.3% | **66.4%** |
| Library 3 (B) | 80.4% | 62.9% | **63.7%** |
| Library 3 (A+B) | 96.5% | 84.5% | **85.7%** |

Movement is under 1.3 points either way, so the headline - nanopore recovers more
designs than Illumina at equal effort, and A+B beats both encodings alone - does
not rest on the depth definition. Both are in `platform_comparison.json`
(`matched_depth_per_block` and `recovered_wells_per_block`).

### 2. Outcome populations replace mean identity

The mean-identity figure is gone, replaced by
`platform_sequence_populations`: the percentage of **consensus sequences**
(nanopore) or **reads** (Illumina) falling in each outcome. Same three
categories and same layout as the recovery figure, so the two read together;
the difference is the denominator - observations here, designs there.

| set | nanopore consensuses | | | Illumina reads | | |
|---|---|---|---|---|---|---|
| | Perfect | Screen. | Other | Perfect | Screen. | Other |
| Library 1 (stuffer) | 93.7% | 2.5% | 3.8% | 79.2% | 16.0% | 4.8% |
| Library 3 (A) | 93.7% | 2.4% | 3.9% | 76.9% | 18.7% | 4.4% |
| Library 3 (B) | 92.2% | 2.7% | 5.2% | 78.5% | 16.2% | 5.4% |
| Library 3 (A+B) | 93.0% | 2.5% | 4.5% | 77.7% | 17.5% | 4.9% |

A nanopore consensus is perfect ~93% of the time; a single Illumina read is
perfect ~78% of the time, with most of the difference landing in Screenable.
Averaging over reads inside a well is what buys that, and it is the same point
the mean-identity figure made (consensus 99.7% against read 95.6%) in units that
match the rest of the set. The mean-identity numbers remain in
`platform_comparison.json` under `read_accuracy`.

Note the asymmetry, unchanged from the previous entry: our `mixed_variants`
grade sits in Other and has no Illumina analogue, since a read is one molecule.

### 3. Widths cut to fit the content

| figure | before | now |
|---|---|---|
| `platform_sequences_per_well` | 89 x 47 mm | **68 x 39 mm** |
| `platform_reference_recovery` | 183 x 58 mm | **114 x 50 mm** |
| `platform_sequence_populations` | 89 x 50 mm | **114 x 50 mm** |

The recovery and population figures are now **two panels**: Perfect on a 0-100%
axis, Screenable and Other on a magnified one. On one axis the small categories
were a smear along the floor and two thirds of a 183 mm canvas was blank. Split,
the width goes to data and the small bars become readable.

**The magnified ceiling is computed from the data.** A hard-coded 5% clipped the
Library 3 (A) Illumina Screenable bar at 5.26% the first time it ran - a clipped
bar still looks like a bar, so nothing about the figure said it was wrong. It now
takes a round ceiling above the tallest bar in the panel.

`platform_read_accuracy.{pdf,png,svg}` was deleted rather than left to rot beside
the figures that replaced it.

### Added

10 more tests in `tests/test_platforms.py` (193 pass), covering the
apportionment (single-block plate, shared plate, exact sums under rounding,
allocated exceeding recovered) and the population counts.

### Next steps

Unchanged from the previous entry, minus the depth question, which is settled:

1. **`SUMO_B` looks weaker than A on two independent measures** - 6.7% empty
   wells against 2.8%, and the lowest Illumina read accuracy at 98.4%. Worth
   understanding before the encodings are treated as interchangeable.
2. **Audit `figures.py::_save`** for the tight-bbox issue; `fig1`-`fig4` still
   inherit it, while `save_figure` does not.

---

## 2026-08-20 (fifth) — Nanopore clone picking against pooled Illumina

Three figures remaking and extending an earlier comparison, over the same
library sets. Analysis only: no stage, threshold or default was touched.

Library sets, all with 11 culture plates on one PCR barcode:

| set | barcode | culture plates | designs |
|---|---|---|---|
| Library 1 (stuffer) | RP08 | `LAB_P1`-`LAB_P11` | 332 |
| Library 3 (A) | RP05 | `SUMO_A_P1`-`P11` | 342 |
| Library 3 (B) | RP05 | `SUMO_B_P1`-`P11` | 342 |
| Library 3 (A+B) | RP05 | both | the same 342, recovered if either encoding was |

A and B are two DNA encodings of one design set: their 342 design keys are
identical, which is what makes the A+B union meaningful rather than a fourth
library.

### The join is exact

Illumina names an encoding with `|` where ours use `_`; parsing the block prefix
off both gives (encoding, block, design key), and the two datasets then agree
**exactly**: 332, 342 and 342 references with the same block numbers and nothing
left over on either side. Checked before any figure was drawn.

### Matched sampling effort

Picking N colonies from a block and taking N reads of that block are the same
investment, so Illumina is subsampled per block to the number of wells the
nanopore run sequenced of that block: **3,000 reads across 51 blocks** (1,009
stuffer, 1,016 A, 975 B), drawn with `random.Random(141142)` in sorted block
order so the sample does not depend on file or dict ordering.

Without this the comparison only says that a deep pool sees more designs than 96
colonies: at full depth Illumina recovers 97-100% of every set.

**Caveat on the definition.** "Wells sequenced of a block" counts wells that
*yielded* a clone of that block, because a well that yielded nothing cannot be
attributed to a block. That is 975-1,016 rather than the 1,045 allocated, so
Illumina is given ~3-7% less effort than the colonies actually picked. Erring
this way understates Illumina slightly; the alternative would require inventing
an attribution for empty wells.

### Allocated wells: 95 per plate, not 96

n = 1,045 per set = 11 plates x 95 wells. H12 is excluded, and the data says so
rather than the layout: **H12 yielded no clone in any of the 33 culture plates**,
and its reads do not behave like a picked well - RP05 H12 has 40 uniquely
assigned reads against 4,069 in A1, with 6,835 `no_match`. It is a control, not
a colony. Exposed as `--wells-per-plate` rather than hard-coded.

### Figure 1: distinct sequences per allocated well

| set | 0 | 1 | 2 | 3+ |
|---|---|---|---|---|
| Library 1 (stuffer) | 3.5% | 91.0% | 5.2% | 0.3% |
| Library 3 (A) | 2.8% | 91.5% | 5.5% | 0.3% |
| Library 3 (B) | 6.7% | 86.8% | 6.1% | 0.4% |

The source figure read approximately 3/93/4, 2/92/5 and 5/90/5 - the same
picture. "Distinct" is by sequence content (sha256), not by name, so two names
resolving to one sequence count once. Chimeric consensuses count: they are
distinct sequence present in the well.

### Figure 2: recovered references

Perfect, as a percentage of the designed library:

| set | nanopore | Illumina, matched depth | Illumina, full depth |
|---|---|---|---|
| Library 1 (stuffer) | 75.0% | 54.5% | 98.2% |
| Library 3 (A) | 85.7% | 67.3% | 97.4% |
| Library 3 (B) | 80.4% | 62.9% | 97.7% |
| Library 3 (A+B) | 96.5% | 84.5% | 100.0% |

Nanopore reproduces the earlier analysis almost exactly (75 / 86 / 80 there,
75.0 / 85.7 / 80.4 here). Illumina at matched depth is a little below the earlier
57 / 68 / 65, consistent with the slightly smaller depth this definition gives
it. Full-depth numbers are in the JSON, not the figure.

**The A+B result is the headline: two encodings of the same design set recover
96.5% of designs perfectly, against 85.7% and 80.4% separately.** Redundant
encoding buys more than deeper sequencing of one encoding would.

Screenable and Other stay under 3.3% everywhere (`platform_reference_recovery.csv`).

**One asymmetry to keep in mind.** Our `mixed_variants` grade - a polyclonal well
- falls into Other, and Illumina has no analogue because a read is one molecule.
Other is therefore not strictly like for like, though it involves 0-5 designs per
set.

### Figure 3: read accuracy

Mean identity to the design:

| set | nanopore read | nanopore consensus | Illumina read |
|---|---|---|---|
| Library 1 (stuffer) | 96.52% | 99.75% | 99.72% |
| Library 3 (A) | 95.49% | 99.66% | 99.54% |
| Library 3 (B) | 95.74% | 99.72% | 99.01% |
| Library 3 (A+B) | 95.63% | 99.69% | 99.28% |

**The nanopore consensus column is not in the original figure and was added
deliberately.** Raw nanopore reads are 3-4 points behind Illumina reads, and a
figure showing only those two would say the pipeline is less accurate than
Illumina. What the pipeline delivers is the consensus, which meets or beats
Illumina read accuracy on all four sets. Omitting it would have been misleading.
It is a separate bar, not blended into the nanopore one.

Notes: the y axis is truncated at 90% because every value exceeds 95% and a
0-100 axis hides the differences the figure exists to show. Illumina identity is
`1 - edit_distance/reference_length` from their per-read table; nanopore identity
is the edlib alignment identity recorded at assignment - different aligners, same
quantity. Accuracy is insensitive to the subsample (full-depth Illumina means are
in the JSON and differ by <0.15 points).

### Two bugs worth recording

1. **Illumina recovery came out at 115.7%.** Reads were being attributed to a
   library by the *well* they came from (`expected_group`), but a read from a
   SUMO_A well can match a B reference, and those references are not in A's
   universe. Membership now follows the reference a read actually matched, which
   is also the scientifically right answer: such a read is evidence about B.
   A percentage above 100 is a gift - it makes a definition error impossible to
   miss. Had the error gone the other way it would have looked like a result.
2. **The A+B set had zero nanopore reads.** The per-read loop `break`s after the
   first matching library set, so reads matched `sumo_a` and never reached the
   union. A read belongs to *every* set that contains it.

Both are covered by tests in `tests/test_platforms.py` (183 tests pass).

### Added

- `src/nanopore3/platforms.py` — joins, matched-depth subsampling, the three
  quantities. No plotting.
- `scripts/compare_platforms.py` — CLI and figures.
- `figures.use_journal_style()` and `figures.save_figure()` — the house style
  (Arial, black furniture, inward ticks, transparent SVG, untrimmed canvas) now
  lives in one place and `compare_replicate_plates.py` uses it too.
- `tests/test_platforms.py` — 19 tests.

Reproduce:

```
python scripts/compare_platforms.py \
    --run-dir runs/260608-AI-DBTL-v4 \
    --illumina ".../20260724_AI_DBTL_ILLIMINA/Analysis/results/LAB_STUFFER_SUMO_CONCORDANT"
```

Outputs in `runs/260608-AI-DBTL-v4/figures/`: `platform_sequences_per_well`,
`platform_reference_recovery`, `platform_read_accuracy` (each PDF, 600-dpi PNG
and transparent SVG at 89 or 183 mm), plus
`platform_reference_recovery.csv` and `platform_comparison.json`.

### Next steps

1. **Ask whether the matched-depth definition should count allocated wells
   instead of recovered ones.** It moves Illumina up slightly; the choice is the
   experiment owner's, and the script would need one extra flag.
2. **`SUMO_B` recovers 6.7% empty wells against A's 2.8%**, and B's Illumina read
   accuracy is the lowest at 99.01%. Both point at the B encoding or its
   handling; worth a look before the encodings are treated as interchangeable.
3. **Audit `figures.py::_save` for the tight-bbox issue** (carried from the
   2026-08-20 (fourth) entry). `fig1`-`fig4` still inherit it; the new helper
   `save_figure` does not.

---

## 2026-08-20 (fourth) — Transparent SVG, and the figure was not 89 mm wide

Output-only change; no number moved.

### Added

`replicate_concordance.svg` alongside the PDF and PNG, saved with
`transparent=True` so it drops onto any background: both the figure and axes
patches are `fill: none`. `svg.fonttype` is set to `none` so text stays text -
editable in Illustrator, and the Arial the figure asks for is the Arial that
renders, rather than glyph outlines. Verified: `font-family: 'Arial',
'Helvetica', 'DejaVu Sans', sans-serif`, 9 kB.

The in-bar text is white but only ever sits on a filled segment, so nothing
disappears against a light background.

### The bug this uncovered

`figures.py` sets `savefig.bbox: "tight"` for every figure in the repo, which
trims the canvas to its content. This figure was therefore being written at
**87.0 mm** wide, not the 89.0 mm single column it was laid out for. A journal
scaling it up to fill the column would have scaled the type with it, which
defeats the entire point of setting 6-7 pt sizes by hand.

Measured after the fix:

| format | width | height |
|---|---|---|
| SVG | 89.00 mm | 39.62 mm |
| PDF (MediaBox) | 89.00 mm | 39.62 mm |
| PNG at 600 dpi | 88.98 mm | 39.62 mm |

Note that `savefig(bbox_inches=None)` does **not** disable it - `None` means
"use the rcParam". The figure sets `savefig.bbox: "standard"` locally instead.

Removing the trim exposed what the trim had been hiding: the `100%` tick label
ran off the canvas and was being silently absorbed. Margins are now explicit
(`subplots_adjust(left=0.16, right=0.955, ...)`).

**Lesson: a tight bounding box makes a figure's physical size an output of its
content rather than a property you set.** For anything going into a fixed column
width, set the margins and leave the canvas alone. The other figures in
`figures.py` are still saved with the tight bbox and are worth auditing for the
same reason before submission.

### Next step

Audit `figures.py`'s `_save` for the same issue if the report figures are to be
submitted at a fixed column width; `fig1`-`fig4` currently inherit the trim.

---

## 2026-08-20 (third) — Figure restyled to the requested house style

Presentation-only change to `replicate_concordance`; no number moved.

| element | before | now |
|---|---|---|
| font | Helvetica first, Arial as fallback | Arial first, embedded as `ArialMT` / `Arial-BoldMT` (checked with `pdffonts`) |
| text | dark grey `#1B1F24` | black, white only where it sits on a filled segment |
| axes and ticks | grey spines, outward ticks | black spines, ticks pointing in |
| bars | white segment edges | black, 0.5 pt, on bars and on every segment |
| in-bar text | left-aligned | centred in the bar |
| bar labels | `Plate 1` over `RP06 \u00b7 SUMO_A_P1` | `Plate 1` alone |
| ratios | `82/82 matched`, `83/84 exact seq.` | `82/82 wells with same content`, `83/84 clones with identical sequence` |
| legend | two rows | one row |

The barcode-to-culture-plate mapping dropped off the bars, so it now lives only
in the caption, `replicate_concordance.json` and the notebook. Anyone reading the
figure alone cannot tell which barcodes "Plate 1" refers to - state it in the
caption.

Legend labels shortened to fit one row: "Same content, identical seq.", "Same
content, different seq.", "Pooled only", "Dedicated only", "Different content".

**One row is measured, not assumed.** The three current labels fill 89 mm at 6 pt
almost exactly, so a fourth category would have run off the page silently. The
figure now estimates the required width from label lengths and drops to two rows
if a single row would need less than 5 pt. Percentages were kept in the in-bar
text, as in the source figure; say so if they should go.

---

## 2026-08-20 (second) — Figure labels name their unit

Label-only change to `replicate_concordance`; no number moved.

The two in-bar ratios count different things and the first draft did not say so,
which was immediately confusing when read:

- `82/82 wells, same clone` — **wells**, over the wells the dedicated barcode
  saw. Tests demultiplexing, assignment and deconvolution: did the pool put the
  same clone in the same well of the same culture plate.
- `83/84 clones, identical seq.` — **clones**, over the clones both datasets
  found. Tests consensus base-calling, and is a finer check applied inside the
  matched set.

The denominators differ because a polyclonal well holds more than one clone: on
plate 1, 82 wells contain 84 comparable sequences (A09 holds two designs, F08 a
design and a chimera); on plate 2, 93 wells contain 96 (A03, B01 and F05 each
hold two designs). Well H04 sits in both rows at once - matched as a well,
not identical as a sequence - which is only readable if the units are named.

Legend labels changed with them: "Matched, exact seq." became "Same clone,
identical seq." because segments are wells, and "matched" alone read as a
sequence match.

**Lesson: a ratio in a figure needs its unit in the label, not in the caption.**
Two ratios stacked in one bar with different denominators and no unit named is a
misreading waiting to happen, and the reader has no way to resolve it from the
figure alone.

---

## 2026-08-20 — Replicate concordance: two culture plates sequenced twice

### Why

Two SUMO A culture plates were sequenced twice in the 260608 run: once on their
own barcodes (RP06, RP07) and once inside the compressed RP05 pool, where 22
culture plates share every well. That is an independent measurement of the same
physical wells, so the agreement between them tests demultiplexing, assignment,
consensus and compressed-PCR deconvolution against data rather than against the
pipeline's own summaries. An equivalent figure exists from the legacy analysis,
so this also compares the two pipelines on the same question.

**This is an analysis, not a pipeline change.** No stage, threshold or default
was touched, so no result in `runs/260608-AI-DBTL-v4/stages/` changes.

### The pairing was tested, not assumed

The user's recollection was that RP06 and RP07 correspond to SUMO A plates 1 and
2. That is an experimental claim, so `replicate.verify_pairing` checks it from the
data: every design carries a block name, and the run's own
`compressed_pcr.blocks` map says which culture plate that block came from.

| barcode | blocks seen | culture plate implied |
|---|---|---|
| RP06 | `A_Block_1` only, 83/83 clones | SUMO_A_P1 |
| RP07 | `A_Block_2` only, 95/95 clones | SUMO_A_P2 |

Unanimous, with nothing unresolved. `compare_replicate_plates.py` refuses to plot
an unconfirmed pairing unless `--allow-unconfirmed-pairing` is passed, so a
mistaken pairing fails loudly instead of producing a plausible figure.

### What was added

- `src/nanopore3/replicate.py` — comparison logic, no plotting: `load_calls`,
  `compare_wells`, `summarise`, `verify_pairing`. Chimeras participate, keyed by
  their parent signature, because a chimera is real sequence in a real well.
- `scripts/compare_replicate_plates.py` — CLI and figure.
- `tests/test_replicate.py` — 13 tests. 164 pass overall.

Reproduce with:

```
python scripts/compare_replicate_plates.py \
    --run-dir runs/260608-AI-DBTL-v4 \
    --pair RP06=RP05/SUMO_A_P1 \
    --pair RP07=RP05/SUMO_A_P2
```

Outputs into `runs/260608-AI-DBTL-v4/figures/`: `replicate_concordance.pdf/.png`,
`replicate_concordance_wells.csv` (one row per well), `replicate_concordance.json`
(summary plus the pairing check).

### Result

| | RP06 / SUMO_A_P1 | RP07 / SUMO_A_P2 |
|---|---|---|
| dedicated wells matched | 82/82 (100%) | 93/93 (100%) |
| shared clones byte-identical | 83/84 (98.8%) | 96/96 (100%) |
| wells recovered only from the pool | 12 | 1 |
| wells recovered only from the dedicated barcode | 0 | 0 |
| clones called as different designs | 0 | 0 |
| chimeric clones matched | 1 | 1 |

Every well the dedicated barcode saw was also recovered from the pool, assigned
to the right culture plate, with the same clone. Legacy reported 82/83 and 93/95
matched on the same plates with one well not called and one pool-only; the new
analysis has no unmatched wells.

**Definitions.** The bar spans every well with a sequenced clone in *either*
dataset, so a well only the pool recovered widens the bar instead of vanishing.
The matched fraction is over wells the *dedicated* barcode saw — a pool-only well
is extra yield, not a disagreement, and putting it in the denominator would have
reported plate 1 as 87% concordant when nothing disagreed.

### Two findings

**1. Row G of RP06 failed at the bench, and the pool rescued it.** All 12
pool-only wells on plate 1 are G01-G12. Reads per row on RP06:

| row | A | B | C | D | E | F | **G** | H |
|---|---|---|---|---|---|---|---|---|
| reads | 4,898 | 3,609 | 3,283 | 3,506 | 3,329 | 3,524 | **37** | 3,152 |

37 reads against ~3,000-4,900 elsewhere, and RP07's row G is normal at 2,947, so
this is that plate's colony PCR or pooling, not the pipeline. The compressed pool
recovered all 12 clones. This is an argument *for* compressed PCR: the redundant
route covered a row the dedicated route lost.

**2. The single non-identical clone is one ambiguity code, not a disagreement.**
Well H04, design `A_Block_1_dTF085_88_0`, 30 reads each side, 240 bases, one
position differs: `N` on RP06 against `G` on RP05. RP06 called that well
`mixed_variants` (polyclonal) while the pooled reads reached the 0.60 support
threshold for G. No base is called differently anywhere in 180 shared clones.

### How strong is this evidence

The consensus is reference-guided, so two runs of a clone that matches its design
will both return the design sequence, and 171 of the 179 exactly-agreeing
designed clones are graded `perfect`. The informative cases are the 8 that are
not: 3 `screenable`, 2 `premature_stop`, 1 `frameshift` and both chimeras agreed
byte-for-byte from independent read sets. Deviations from the design reproduce,
which they would not if consensus were smoothing reads toward the reference.

The claim this supports strongly is that demultiplexing, assignment and
deconvolution are right. The claim it supports weakly is per-base accuracy,
since most clones are perfect matches where agreement is close to trivial.

### Figure

`figures/replicate_concordance.pdf`, 89 mm single column, print style from
`figures.use_print_style`.

The palette needed changing rather than copying. The source figure separates
"matched, different seq." (pale teal) from "not called" (pale grey), and

```
scripts/validate_palette.py "#2E6F6A,#93C4BD,#B4462F,#D9A441,#DCE3E7" --pairs all
```

puts that pair at dE 13.6, below the 15.0 normal-vision floor — those two thin
slivers genuinely are hard to tell apart, which matters when a sliver is the
entire finding. The chosen set `#2E6F6A, #6FA9A0, #A8442E, #E0A33A` (plus the
neutral `#EEF1F3`) reaches worst normal-vision dE 17.8 and worst CVD dE 9.5,
above the 8.0 target. The validator's per-colour chroma and contrast checks are
relaxed here on purpose: these are large area fills carrying a legend and white
in-bar text, not thin line series, and neutrals are already exempt by the
convention in `figures.py`.

Legend labels are clearer than the source's rather than verbatim: "Pooled
barcode only" for what legacy called "Polyclonal Seq only", and "Dedicated
barcode only" for "Not called", both of which are ambiguous out of context. Say
so if the verbatim wording is needed for a direct comparison.

### Next steps

1. **RP08 has no dedicated replicate**, so the 11 LAB culture plates have no
   equivalent check. If any LAB plate was sequenced separately, add it as a
   `--pair` — the script takes any number.
2. **Ask at the bench what happened to row G of the RP06 plate.** The pipeline
   has nothing more to say about it; 37 reads is an absence of data.
3. **Consider running this comparison as a pipeline stage** when a config
   declares a dedicated replicate of a pooled culture plate. It is currently a
   script because the pairing is knowledge about the experiment, not about the
   data.

---

## 2026-08-19 (fifth) — QC output follows the culture plates, and chimeras are graded

### What changed

Three changes, all **scientific** in the sense that they change what appears in
the delivered outputs, though none change how any read is assigned or any
consensus base is called.

1. **The consensus tree nests by culture plate where the layout resolves it.**
   Previously every plate barcode wrote `<plate>/<well>/`. When
   `compressed_pcr` is configured, a clone's block names its source culture
   plate, so RP05 now writes `RP05/SUMO_A_P6/A01/...` and RP08 writes
   `RP08/LAB_P3/B07/...`. The deconvolution was already being computed and
   reported; it just was not reaching the thing a person actually opens.
   Barcodes with no pooling (RP01-RP04, RP06, RP07, RP09-RP11) stay flat at
   `<plate>/<well>/`, and so does any clone whose culture plate does not
   resolve to exactly one source — an unresolved set is joined with `|` and
   deliberately does **not** become a directory, because a directory named for
   a set of plates would read as a claim the data does not support.

2. **Chimeric clones are graded and written into the QC tree.** They are real
   sequence present in a real well, so excluding them from QC was hiding data.
   `04b_chimera` now emits `scaffolds.fasta` keyed by a `chimera_id`
   (`chim-<digest>`) and resolves each clone's culture plate from its first
   parent's block. `05_qc` grades each chimeric consensus against its own
   spliced scaffold — not against any single parent, which would score a
   correct chimera as badly truncated — and appends the rows to `qc.csv.gz`.
   Chimera files are marked in three places so they can never be mistaken for
   a designed clone: `chimera__` in the filename, `kind=chimera` in the FASTA
   header, and a `kind` column in `index.csv`. PCR-origin chimeras remain
   unwritten, per the 2026-08-19 (third) entry.

3. **The HTML report renders structured sections.** See below.

### The report bug, and why it is worth an entry

Run v4 completed stages 01-05 and then died in `06_report`:

```
TypeError: int() argument must be a string, a bytes-like object or a real number, not 'dict'
```

`write_html_report` took `Mapping[str, Mapping[str, int]]` and did `int(value)`
on every cell. The culture-plate section I had added is a mapping *of mappings*
and contains a string (`weakest_culture_plate: SUMO_B_P11`). The contract was
implicit — no type check runs on this path — so adding a legitimate section
turned into a crash in the final stage of a twenty-minute run.

The fix is not to flatten my data to fit the renderer. `_rows` now recurses one
level for nested mappings and `_format_value` handles bools, ints (thousands
separators), floats and strings, escaping all of them. Four tests in
`tests/test_science.py::ReportRenderingTests` cover flat counts, a nested
per-plate section, a string value, and HTML escaping.

**Lesson: a renderer that coerces its input is a landmine for whoever adds the
next section.** The type annotation said `int`, nothing enforced it, and the
cost landed at the end of the longest stage rather than at the call site. When
a function's declared input type is not checked, prefer widening the function
over narrowing the caller.

### Run 260608-AI-DBTL-v4

Resumed after the fix; stages 01-05 were reused unchanged (config digest
`e5d6c132...`) and only `06_report` re-ran. Outputs in
`runs/260608-AI-DBTL-v4/`.

Consensus tree, from `stages/05_qc/consensus_by_plate/index.csv`:

| | files written |
|---|---|
| nested `<plate>/<culture>/<well>/` | 3,194 (RP05 2,126 · RP08 1,068) |
| flat `<plate>/<well>/` | 5,649 |
| total | 8,843 |

RP05 recovers all 22 pooled culture plates, RP08 all 11 — confirmed by
directory listing, not by the summary that computed it.

Chimeric consensuses now in QC: **310** of 21,094 graded entries. Their grade
distribution is nothing like the designed clones', which is the expected
result and a useful sanity check:

| grade | chimeric | designed |
|---|---|---|
| frameshift | 182 | 176 |
| truncated | 103 | 50 |
| premature_stop | 11 | 35 |
| perfect | 3 | 7,759 |

A chimera is a junction between two designs, so a frameshift or a truncated
alignment against its own scaffold is the norm; 3 perfect ones are cases where
the splice happens to be in-frame and clean. If chimeras had come out looking
like designed clones, the scaffold synthesis would have been suspect.

### Next steps

1. **RP03/RP04 still produce a disproportionate number of chimeras** — 187
   clones from 30.6% of reads against RP05's 68 from 35.5%, on the same
   opTF001 library. Unexplained. RP04 alone wrote 4,162 consensus files, far
   more than any monoclonal plate, which is consistent with polyclonal wells
   but has not been checked against the intended layout (unknown at the bench).
2. **SUMO_B_P11 built 42 sequences against a median of 96** across the other
   21 culture plates, flagged in `fig4_culture_plates`. This is a bench
   question — colony density, or a plate that was under-picked.
3. **`low_depth` is 12,251 of 21,094 groups (58%).** This is the dominant
   yield loss and is a reads-per-well problem, not a pipeline one. Worth
   deciding whether `minimum_depth: 6` is the right threshold for this depth
   of sequencing, and reporting the yield curve if it is changed.
4. **`motif_missing` runs 12-17% on RP05-RP11 against ~5% on RP01-RP04.** The
   split is by construct, so the per-library `reverse_motif` overrides are the
   first place to look.
5. **Lint debt** (carried from earlier entries): the over-long CSS line in
   `report.py` is fixed, but `ruff check src/` still reports pre-existing
   findings and CI does not gate on it.

---

## 2026-08-19 (fourth) — Culture-plate recovery, and a metric that was 2.9x wrong

### What was added

Deconvolution was only visible as a read-fate bar, which says whether a read
resolved but not *which* culture plate it came from — the point of compressing
plates in the first place. `fig4_culture_plates` adds, per pooled colony-PCR
plate, a 96-well map of how many distinct source plates each well recovered and a
bar per source plate. The report gains a **Culture plate recovery** section.

The pooled denominator comes from `compressed_pcr.pcr_plates`, not from what was
observed: a culture plate contributing nothing must show as missing rather than
silently shrinking the denominator from 22 to 21.

### The metric was wrong, and the experiment owner caught it

The first version reported 200-375 "clones" per culture plate. The experiment
owner questioned it immediately: RP05 has 96 wells and each well draws one colony
per culture plate, so **96 is the ceiling** — 375 is impossible.

The cause: it counted consensus **groups**, including 3,769 that fell below the
depth floor and produced no sequence at all. Counting groups inflated recovery by
**2.9x** and put more sequences on a culture plate than the plate has wells.

Corrected to count consensus sequences actually built:

| | RP05 | RP08 |
| --- | --- | --- |
| consensus sequences | 2,058 | 1,052 |
| groups below depth (not credited) | 3,769 | 2,321 |
| per culture plate, median | **97** | **96** |
| source plates with a consensus per well | median 21 of 22 | median 11 of 11 |

**97 and 96 against a ceiling of 96** is what near-complete recovery looks like:
essentially every (well, culture plate) pair yields one consensus. The occasional
102 or 106 is a well where one culture plate contributed two distinct designs.

The per-well figure also became more honest: median 21 of 22, not 22 of 22,
because a source plate is now only credited when it produced a usable sequence.

**SUMO B plate 11 built 42 sequences against a median of 96** across the other 21
— under half. Worth checking at the bench.

### Two figure defects found by rendering

The outlier highlight compared against the **mean**, which the outlier itself
drags down: at 42 against a mean of 94 it failed the 0.4x test and went
unhighlighted. Comparing against the median catches it. Titles also collided once
they grew, and were shortened.

### Lessons

**A count needs a ceiling check.** "200-375 per culture plate" was reported
without asking what the maximum possible value was. One division — 96 wells, one
colony each — would have caught it before it reached a figure.

**"Groups" and "sequences built" are not interchangeable**, and on this run they
differ by 2.9x because most groups are too shallow. Any per-plate or per-well
count has to say which it means.

**An outlier detector must not use a statistic the outlier distorts.** Flagging
against the mean hid exactly the case the flag exists for.

### Next steps

1. Check SUMO B plate 11 at the bench: 42 sequences against a median of 96.
2. Consider showing groups-below-depth per culture plate too, since 3,769 of
   5,827 RP05 groups are shallow and that ratio may itself vary by plate.
3. Unchanged from earlier: investigate RP03 and RP04's chimera rate.

---

## 2026-08-19 (third) — Chimera detection wired in and run across all plates

Run `runs/260608-AI-DBTL-v3`, **20.8 minutes** including the new stage.

**`04b_chimera` is now a pipeline stage**, beside consensus rather than inside it
so a chimera setting change does not invalidate the reference-guided
consensuses, and named `04b` so existing stage numbering — and the QC, report,
export and figure code — stays untouched.

**PCR-origin chimeras are counted but not written.** `classify_origin` moved from
the script into `nanopore3.chimera`, because it now decides *what gets written*
rather than only how output is labelled, and that belongs in tested code.
`write_pcr_origin` restores them for inspection. New `chimera` config section,
disabled by default so other profiles are unaffected.

### Results across all 11 plates

| | |
| --- | --- |
| chimeric reads | 36,419 of 1,174,553 with an insert |
| **clones written** | **310** (all assembly-origin, all two-parent) |
| PCR-origin groups counted, not written | 74 |
| declined to scaffold (three parents) | 3 |

Median 33 reads per clone, up to 579, and **median 0 ambiguous bases** — the
chimeric consensuses are clean, which is what a real clone should give.

| plate | clones | share of assigned reads |
| --- | --- | --- |
| RP04 | **187** | 30.6% |
| RP05 | 68 | 35.5% |
| RP03 | 23 | 2.1% |
| RP08 | 16 | 15.6% |
| RP09 | 11 | 2.2% |
| RP01, RP10, RP11 | 0 | 2.0-3.5% |

**RP04 and RP03 are out of proportion to depth.** RP04 gives 187 clones from
30.6% of reads while RP05 gives 68 from a larger 35.5%; RP03 gives 23 from 2.1%.
Both are opTF001, so this is a library or plate effect rather than a depth
effect, and is worth investigating at the bench.

### Agreement with the legacy analysis on RP05

Of legacy's 28 chimeras, **25 are written**, and the three that are not are each
correct behaviour rather than a miss:

- **E10 and F02** have parents in different assembly blocks, so they are PCR
  template-switch artefacts and are deliberately withheld. Legacy wrote them as
  clones.
- **D03** has three parents, where `synthesise_reference` declines to splice
  because the junctions are not independently located well enough. A genuine gap.

The stage reproduced the standalone script's RP05 figure exactly (68
assembly-origin) across a separate full run.

### Lessons

**Where a decision lives matters as much as whether it is right.** Origin
classification was fine in a script while it only labelled output; the moment it
decided what to write, leaving it there would have put a scientific filter
outside the test suite.

**A disagreement worth having is one you can explain per case.** Three legacy
calls are missing and each has a specific reason — two artefacts withheld by
design, one real limitation. More useful than a matching count.

**A failed `&&` link silently skipped a lab-notebook entry.** The previous
session's chimera entry was never written, because a lint failure short-circuited
the chain while the `git commit` on the following line still ran. The record and
the commit disagreed. Notebook writes should not ride on a chain that can fail.

### Next steps

1. **Investigate RP03 and RP04's chimera rate**, far out of proportion to their
   read share and specific to opTF001.
2. Recover the three-parent case: locate multi-junction breakpoints well enough
   to splice, or make the reference-free fallback produce a usable consensus.
3. Grade chimeric consensuses as the reference-guided ones are graded, so
   `perfect`/`screenable` applies to them too.
4. Surface chimera counts in the HTML report and as a figure.

---

## 2026-08-19 (second) — Chimeric clones are recovered instead of discarded

### Validation against the legacy analysis

**Sequence agreement is essentially exact.** 2,057 of our 2,058 RP05 consensuses
appear in the legacy set (Jaccard 94.9%), 1,865 graded perfect on both sides.
**97.4% of legacy consensuses are exactly ours plus `AGTGGA`** — the vector
overhang — because the legacy consensus is not trimmed to the reference span
while ours is. After removing those 6 bases, edit-distance-to-reference agrees on
97.6%. The 48 cases where legacy said screenable/other and we said perfect are all
explained by it: legacy checks 6 constant bases we do not.

**The 110 consensuses only legacy has are entirely a depth choice**: median 1
read, maximum 5, all under our `minimum_depth: 6`.

### A claim I made and then disproved

I wrote that the 28 wells legacy flagged as chimeric were "reported as ordinary
assignments by us". Checking showed otherwise: of the 57 parent designs legacy
named, **46 never reached consensus at all**, 4 were low_depth, 4
`mixed_variants`, 3 `consensus_pass`. The chimeric material was being *discarded*,
not mis-reported. Better than I said, but still a loss.

`mixed_variants` is **not** a chimera call. It says the reads assigned to one
design disagree with each other beyond the support threshold. A chimera can cause
that, but a well holding *only* chimeric molecules yields a clean consensus that
is confidently wrong, and `mixed_variants` has other causes entirely.

### What was built

The experiment owner's proposal: match k-mers positionally, then group reads
sharing a pattern and build a consensus per group. Better than the read-to-read
clustering in `nanopore3.clustering`, because it compares read **segments to
references** rather than reads to each other, avoiding the regime measured on
2026-08-13 where clones 94% identical cannot be separated at nanopore error rates.

`src/nanopore3/chimera.py` walks each insert in 90 nt windows at 45 nt steps,
matches every window independently, and collapses the calls into a signature. A
clean clone gives one segment; a chimera gives a run of each parent, naming both
and locating the junction. Reads sharing a signature are the same clone, so the
usual depth floor applies: a signature carried by fewer than `minimum_depth` reads
is one noisy window, not a clone.

### The distinction that matters

On RP05, 126 chimeric clones were recovered from previously discarded reads,
including **27 of legacy's 28 exactly**. They split by whether the parents share
an assembly block:

| origin | clones | median reads | in legacy | new |
| --- | --- | --- | --- | --- |
| assembly (same block) | 68 | **140** | 25 | 43 |
| pcr (different blocks) | 58 | **8** | 2 | 56 |

Two designs of one block reach a colony-PCR well only by being the same colony,
so a same-block chimera is a genuine assembled clone. Designs from different
blocks arrive as separate colonies pooled into the shared well, so a chimera
between them can only have formed in the tube. **The 17x difference in read depth
confirms it independently**: real clones are clonally amplified, template-switch
products appear late and stay rare.

---

## 2026-08-19 — 260608 AI_DBTL run: six libraries, eleven barcodes, 2.55 M reads

Run `runs/260608-AI-DBTL-v2`, **13.3 minutes** for 2,546,670 reads.

### Setting up

`/Volumes/TJAF/260608_SUMO_LAB_BS_FS_Ph/AI_DBTL.fastq` (8.4 GB) against six
designed libraries. Everything below was **derived from the reads rather than
assumed**, after two earlier sessions where assuming cost real time.

Two reference files each held **two** sub-libraries under `A|`/`B|` prefixes with
independent block numbering (SUMO, EBM); they were split. EBM references were
"flanked" and were stripped to bare inserts. Library routing per plate barcode was
confirmed empirically at 65-94% before any config was written.

Corrections to what was initially believed:

- **RP04 is opTF001, not SUMO B.** Assumed from abundance; the data said otherwise.
- **RP05 carries both SUMO A *and* B** — 13,124 and 15,223 reads, all 18 blocks
  each, separated at identity margin 0.176. Two complete sets, 22 culture plates.
- **RP06, RP07, RP01-RP03 are real**, not the background I first dismissed them
  as at ~2% of reads; each assigns to its library at 74-97%.
- **RP09, RP10, RP11 are real too** (sensors, EBM, dTF141/142). An earlier entry
  called ~6% of reads "carryover" on the assumption they were unused. **Only RP12
  is unused, and it produced zero reads**, so the plate-barcode false-positive
  rate is effectively nil.

### Settings, and how they were chosen

`max_edits: 6` was swept 2-8 against the RP12 negative control. The
REAL/ABSENT ratio is flat across the range (11.5 to 11.9), so tightening buys no
specificity and only loses yield. That confirms the previously inherited value on
*this* run's own control rather than by inheritance.

### One code change was required

Libraries built on **different constructs** cannot share one motif pair: sensors
and ebm_ab end at a different constant region, so a single global reverse motif
would extract the wrong span. `ReferenceSettings` gained optional per-library
`forward_motif`, `reverse_motif`, `motif_max_edits`, `qc_upstream_constant` and
`qc_downstream_constant`, each falling back to the run-level value.

Either motif may be overridden alone. Requiring both together — as the first
version did — would force a library that only *ends* differently to restate the
shared 5' motif, inviting the two copies to drift apart. Also worth recording:
the supplied 5' motif for sensors/EBM turned out to be literally `CTT` + the
standard one, so no forward override was needed at all.

### The failure, and the fix that matters

The first attempt **crashed 25 minutes in, at QC, with `KeyError: 'A'`**. The
`A|`/`B|` prefix collided with `|`, which the pipeline reserves as its
alias-group separator: `reference_ids` is a `|`-joined list that is split back to
recover the group, so `A|Block_1_...` parsed as alias `"A"`.

A bad reference ID survived demultiplexing, assignment **and** consensus before
failing on re-parse. That is exactly the late, expensive failure this pipeline is
meant to prevent, so `read_reference_libraries` now **rejects `|` in a reference
ID at load**, naming the reason. It fails in preflight, in seconds.

Two process notes from the same episode:

- Patching an already-patched YAML file mangled the block keys twice
  (`"A|"_Block_1"`). Regenerating from the generator worked first time.
- The resume **would have been correctly refused**: the config digest changed
  (`4adb0a82` to `454c16bf`), so the immutability contract would have rejected
  reusing the assignment stage built on the bad names. Verified rather than
  assumed. A fresh run ID was used so provenance matches the final config exactly.

### Results

| Stage | wall | outcome |
| --- | --- | --- |
| 02_demux | 9.9 min | 1,419,291 assigned (55.7%) |
| 03_assignment | 2.7 min | 1,074,747 `assigned_unique` (75.7%) |
| 04_consensus | 0.6 min | 8,269 pass, 12,251 low_depth, 264 mixed_variants |
| 05_qc | <1 min | 8,153 pass (95.5% of evaluable), 380 fail |

Graded tree: **7,759 perfect**, 242 screenable, 264 mixed_variants, 176
frameshift, 50 truncated, 35 premature_stop, 7 mismatched, 12,251 low_depth.
8,533 FASTAs written across 11 plate directories.

**Deconvolution: zero `unexpected_block` across all 1.42 M reads.** All 22 RP05
culture plates and all 11 RP08 plates recovered. SUMO A/B split 46.2% / 53.8%,
consistent with two complete sets pooled.

Per-plate `assigned_unique` is uniform at 71-85% across all eleven barcodes.

### Storage, not CPU, was the earlier bottleneck

Demultiplexing ran at **4,324 reads/s against 786 on the 20260506 run** — 5.5x
faster. Same code, same worker count. That run read its input from OneDrive; this
one from a local external drive. It confirms the I/O diagnosis and settles it: a
throughput figure is meaningless without saying where the input lived.

**Demux was byte-identical between the failed and successful runs** (SHA-256
`4c6d479c`), confirming determinism across 8 process workers on 2.55 M reads.

### Lessons

**Derive the layout, then ask only what cannot be derived.** Every belief that was
wrong here — RP04's library, RP05 holding both halves, RP09-11 being real — was
correctable from the reads in minutes. Three sessions of assuming preceded this.

**A reserved character needs a guard, not a convention.** `|` was load-bearing in
four output paths and documented nowhere enforceable. The convention held until
someone put it in a name.

**Regenerate, do not re-patch.** Two failed `sed`/regex passes over an already
edited config produced malformed YAML both times.

### Next steps

1. **`motif_missing` is 12-17% on RP05-RP11 but ~5% on RP01-RP04.** RP01-04 are
   opTF001 and the rest span SUMO/LAB/sensors/EBM, so this looks
   construct-specific rather than random. Investigate before quoting yield.
2. **`low_depth` is 58.9%**, against 35% on the 20260506 run. Expected given RP05
   compresses 22 culture plates into 96 wells, but it means most well/design
   combinations have under six reads. Consider whether `minimum_depth: 6` suits
   this level of compression.
3. `plate_ambiguous` is 6.2% against 0.1% previously. The direction is expected
   with 11 barcodes rather than 7, but the magnitude deserves a check.
4. RP03/RP04 culture-plate layouts are still undeclared and report at well level.
5. Run `scripts/detect_assembly_errors.py` over this run's demultiplexed reads
   for a chimera and missing-fragment rate on these libraries.

---

## 2026-08-19 — CORRECTION: RWV1-RWV4 were already registered

### What was wrong

For the 20260818 GG HiFi run I registered a new barcode family
`tf544_547_core` holding the bare 24 nt variable cores of the supplied TF544-TF547
primers. That was redundant and worse than what already existed.

The experiment owner pointed out that these runs use **RWV1-RWV4**, already in
`configs/barcodes/reverse_primer_families.csv` as family
`rwv_pet_whole_vector`. Checked: the registered RWV sequences are the **reverse
complements** of the supplied TF primers, differing only by a 1-2 nt G-run length
in the constant tail (`AGCGGGGGATACGGTT` against `AGCGGGGATACGGTT`).

| Supplied | Registered | edit distance over ~59 nt |
| --- | --- | --- |
| TF544 | RWV1 | 1 |
| TF545 | RWV2 | 2 |
| TF546 | RWV3 | 1 |
| TF547 | RWV4 | 2 |

`tf544_547_core` is removed and the profile now uses `family_id:
rwv_pet_whole_vector` with the previously tuned plate settings (trim 12, window
400, max_edits 6, `legacy_unique_threshold`).

### Why the registered family is better, beyond not duplicating

My subset family held only the four barcodes in use, on the reasoning that a read
could not then be mis-called as an unused barcode. That reasoning cost a
**negative control**. With the full eight-member family:

| Barcode | Primer | Plate / round | share of gated reads |
| --- | --- | --- | --- |
| RWV3 | TF546 | plate 2, round 2 | 45.2% |
| RWV2 | TF545 | plate 2, round 1 | 27.1% |
| RWV4 | TF547 | plate 1, round 2 | 9.7% |
| RWV1 | TF544 | plate 1, round 1 | 3.0% |
| RWV5-RWV8 | not used | — | **0%** |

**RWV5-RWV8 receive exactly zero calls.** That is direct evidence the panel and
the experiment agree, and it is unobtainable from a panel containing only the
barcodes expected to appear. Registering the bare cores also discarded the primer
flanks, changing what `trim_bases` means and losing the context that made the
tuned settings transferable.

Also visible: **plate 1 is heavily under-represented**, 12.7% against plate 2's
72.3%, and RWV1 at 3.0% is weak. Worth checking at the bench.

### Lessons

**Search the registry before extending it.** The barcodes were already there under
a different name because they were registered in the orientation the vector uses,
not the orientation the primer order sheet lists. Matching on the variable core
alone, in both orientations, would have found them in one query.

**A panel restricted to what you expect cannot surprise you.** Keeping the unused
members is what turns a barcode panel into a control.

### Status of the 20260818 run

Still **BLOCKED on references**, unchanged by this correction: 98.1% of gated
reads share under 1% of their 15-mers with the 62 supplied plasmids. Barcodes,
read structure and the length gate are all sound; the reference set does not
correspond to this flow cell.

---

## 2026-08-13 (third) — Whole-vector amplicons, and reference-free clustering

Two questions from the experiment owner: would the pipeline handle a 6 kb
whole-vector amplicon with a 300-1000 nt variable insert, and could reads in a
well be clustered without a reference.

### Whole-vector amplicons: it runs, but the thresholds stop meaning anything

Measured on synthetic 200-design libraries of whole vectors.

**Speed is fine.** 4 ms per read, 264 reads/s against 2,955 reads/s for the 1.2 kb
amplicons: about 11x slower per read, so roughly 1.6 hours for 1.5 M reads. The
specificity-weighted k-mer index handles the shared backbone exactly as intended,
because backbone k-mers are owned by every reference and are dropped by
`max_kmer_owners`, leaving only insert k-mers carrying signal.

**The identity threshold becomes meaningless.** A *perfect* read of the **wrong**
design scores **0.9718** identity against a 5.7 kb vector with a 300 nt insert,
far above the `minimum_identity: 0.80` floor. With 95% of the amplicon invariant,
identity measures the backbone, not the design.

Everything therefore rests on the identity margin, which scales with the variable
fraction:

| Insert | read error | worst margin | headroom over the 0.02 threshold |
| --- | --- | --- | --- |
| 300 nt in 5.7 kb | 6% | 0.0242 | 21% |
| 300 nt in 5.7 kb | 10% | **0.0207** | **3.5%** |
| 1000 nt in 6.4 kb | 10% | 0.0664 | comfortable |

No wrong calls occurred, but the synthetic inserts were random and therefore
maximally distinct. A real library of related designs would shrink the margin
further. **The principled fix is to score identity over the variable region
rather than the whole amplicon**, which is not yet implemented.

Reporting accuracy across the entire read length already works: QC aligns the
whole consensus globally, so widening the length gate and supplying whole-vector
references is enough for that part.

### Reference-free clustering

New `src/nanopore3/clustering.py`. Reads in a well are projected onto a seed,
positions where reads systematically disagree are found, and each read is reduced
to its alleles at those positions alone. The invariant backbone contributes
nothing, which is what makes it insensitive to how much sequence the clones share.

Three findings, each from a failure:

**1. The seed must be polished first.** With a raw read as scaffold, its own 5-10%
errors looked like variant positions to every other read: variants were scattered
from position 6 to 2691 when the insert sat at 1200-1500. One round of majority
polishing moved **98% of variant positions inside the insert** and dropped
within-clone distance from 0.270 to 0.233.

**2. Distances need a comparable-position floor.** An indel near a variant leaves
it uncalled, and two reads overlapping at a handful of positions give a mismatch
ratio driven by chance, which widened both tails until they overlapped.

**3. A single greedy pass mis-assigns.** An early centre absorbs reads from
another clone. Added k-means-style refinement: recompute centres from members,
reassign, iterate to a fixed point.

Also: a split is only accepted when **two** clusters reach the size floor.
One cluster reaching it is one clone whose reads scattered, and splitting there
silently discards the shortfall, which reads as low depth rather than a failed split.

Verified: a monoclonal well of 20 reads yields **one cluster of 20 whose consensus
is 100.00% identical** to the truth; three clearly distinct clones yield three
pure clusters.

### The measured limit, and a claim I had to withdraw

I wrote that clustering failures would be conservative — an under-split flagged
downstream as `mixed_variants`. **A test disproved that**: a cluster formed
containing reads from two different clones. That is mis-assignment, not
under-splitting.

The regime is now measured. Three designs sharing 2.4 kb of backbone and
differing across a 300 nt insert are **94.15% identical**, so their true
divergence (6%) *equals* the read error rate (6%). The within-clone p90 and
between-clone p10 distances both land on 0.495 — the distributions touch, and no
threshold separates them. Clones differing across a 900 nt insert separate cleanly.

This matters directly for the whole-vector question, where the variable region is
5-15% of the amplicon: **that is the hard regime, not the easy one.** The limit is
asserted as a test so it cannot regress silently, and documented at the top of the
module. What would fix it is selecting variant positions by co-variation rather
than by frequency alone, which is proper haplotype phasing.

### Lessons

**A convenient claim about a failure mode needs the same test as a feature.**
"Failures are conservative" was plausible, and wrong, and would have been
believed. One test settled it.

**Do not tune a test green.** The mixing test could have been made to pass with a
kinder seed. Encoding the limit as its own assertion keeps the weakness visible.

### Next steps

1. Region-aware identity scoring for whole-vector references: score the variable
   region separately from the backbone, and threshold on that.
2. Co-variation-based variant selection in clustering, to push the separable
   regime below 6% divergence.
3. Wire clustering into the pipeline as an optional reference-free path; it is
   currently a library and is not called by any stage.

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
