# Demultiplex optimisation audit

This audit separates three things that were interleaved in the Nanopore2
notebooks: assay parameters, barcode decision semantics, and execution speed.

## Recovered assay settings

The RP amplicon sweep evaluated 819 combinations on the first 80,000 qualifying
reads from an earlier GEN2 run. It treated RP01–RP04 as expected and RP05–RP12
as off-target; those are plate-level proxy labels, not per-read truth. The
conservative result was `trim=12`, `window=400`, `max_edits=6`: 61,836 expected,
zero off-target, 99 ambiguous and 18,065 unassigned. An F2 rescore preferred
`10/400/10` (64,042 expected, 48 off-target), while F0.5 preferred `10/400/9`
(64,010 expected, 12 off-target). There is no universal winner independent of
the desired precision/yield trade-off.

The well sweep used 37,326 reads from RWV3 in a different run. Columns 1–9 were
treated as expected and columns 10–12 as negative controls, again without
per-read labels. Its F-score rescore selected `12/400/7`: 26,679 expected, 44
control calls, 13 ambiguous and 10,590 unassigned. The original ratio ranking's
`2/50/10` result assigned only 525 reads and is not a useful production optimum.

The current main notebook displays `8/400/6` plate and `10/400/4` well settings
for a later RBNPPh phagemid assay. Its 1,500–3,300 nt length filter is also
inappropriate for the 20260506 FASTQ, whose median is 1,232 nt and q95 is 1,395
nt. Settings are therefore versioned by assay profile, never copied from the
notebook's current kernel state.

## Decision semantics

`legacy_unique_threshold` reproduces the plate/sweep rule: more than one barcode
identity at or below the edit threshold is ambiguous, even if one scores better.
`legacy_unique_best` reproduces the active well code: only a tie among best
barcode identities is ambiguous. `best_margin` is the stricter Nanopore3 mode
with end/orientation conflict evidence and a configured best-versus-second
margin. The 20260506 profile uses threshold-unique for RP and unique-best for
wells; every per-read row retains the selected policy's evidence.

## Parallel execution

Nanopore2 production demultiplexing really was chunked and parallel. It split
FASTQ files into 10,000-read disk shards and selected up to
`min(cpu_count - 2, 16, tasks)` workers. The saved 16-logical-CPU Mac run used 14
threads—not 14 processes—and plate demultiplexing reported roughly 2,500
reads/s. The well stage also used 14 threads and reported roughly 971 reads/s.
The parameter sweeps' similarly named “chunks” were processed sequentially; they
were dataset partitions, not multiprocessing.

That design bounded memory and avoided shared writers, but it also copied large
inputs, created thousands of part files, opened many handles and merged serially.
When output lives in cloud-synced storage, this can become I/O-bound. Nanopore3
instead streams in-memory batches, compiles barcode panels once, bounds pending
work, writes deterministically in the parent and publishes stages atomically.

On the deterministic 20260506 20k pilot, the initial per-read Nanopore3 engine
took 299 seconds for the complete demux stage. After compiled panels, edit
cutoffs, A/C/G/T fast paths and coarse batching, a direct full matcher benchmark
processed 20,000 reads in 50.116 seconds (399.1 reads/s) with two 10k thread
batches. On a 2,000-read worker sweep, serial, four-thread and 14-thread rates
were 370.9, 413.4 and 388.5 reads/s respectively, with identical call counts.
This machine's current thread knee is therefore about 2–4 workers, not all 14.

Use `scripts/benchmark_demux.py` to test a fixed profile and input subset on each
machine. Test serial and threads first. The explicit process backend is intended
for CLI use and uses spawn semantics for macOS, Windows and Linux; it must pass
the platform CI matrix before becoming an automatic default. Avoid nested pools
and enforce `jobs × threads_per_job` within the CPU budget.
