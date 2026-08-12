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
   indexes (~3 libraries, 1,175 references). Check the `empty_read` count in the
   demultiplexing summary before quoting any yield.
2. **Investigate the 973 `motif_missing` and 1,901 `no_match` reads.** That is
   19% of demultiplexed reads discarded. Determine whether they are genuine
   off-target/chimeric molecules or a motif/threshold artefact before treating
   any yield number as final.
3. **Consensus is the next optimisation target** now that it is no longer
   dominated: 1,731 of 2,483 groups fail as `low_depth` and only 286 pass QC.
   That is a scientific yield question (reads per well) more than a speed one.
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
