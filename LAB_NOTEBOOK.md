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

**This change is deliberately not applied.** It depends on a fact only the
experiment owner can confirm: whether `CAGCTT` is a constant linker between the
ATG and every designed insert — in which case trimming it is correct and lossless
— or whether it varies for some designs, in which case trimming would corrupt
them. Evidence is preserved in `runs/20260506-diag-motif-cagctt` and
`runs/20260506-diag-motif-cagctt-e3`.

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
