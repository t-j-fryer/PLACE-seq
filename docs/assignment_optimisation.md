# Reference assignment optimisation

This note records how the `03_assignment` stage works, why it was restructured on
2026-08-12, and what evidence supports the claim that the faster path is
scientifically equivalent. See the [lab notebook](../LAB_NOTEBOOK.md) for the
chronological account.

## The two-tier design

Assignment never aligns a read against a whole reference library. It uses a
cheap k-mer prefilter to propose a handful of candidates, then spends alignment
effort only on those.

1. **Extract the insert once.** The read is reduced to the sequence between the
   configured forward and reverse motifs, in whichever orientation fits. This
   happens once per read, not once per k-mer size.
2. **Shortlist by specificity-weighted k-mers.** `ReferenceIndex` maps each
   canonical k-mer to the alias groups that contain it. A k-mer shared by *n*
   references contributes `1/n`, so a k-mer unique to one reference counts far
   more than a conserved one. `max_kmer_owners` drops k-mers so common that they
   carry no information.
3. **Validate by alignment.** Only shortlisted candidates are aligned with edlib,
   and the identity, coverage, and margin thresholds decide the call.

`kmer_sizes` declares a **cascade**, listed most specific first (`[15, 11, 9]`).
A smaller k is consulted only when the larger ones proposed nothing acceptable,
and only when it proposes a reference the larger ones did not — alignment
evidence does not depend on k, so re-ranking the same candidate set cannot change
the outcome. Alignment geometry is cached per read, so each candidate reference
is aligned at most once regardless of how many k-mer sizes are consulted.

## The rescue policy

If no k-mer size yields an acceptable call, `rescue_policy` decides how hard to
look:

| Policy | Behaviour |
| --- | --- |
| `kmer` (default) | Re-shortlist at the smallest k with no score floor, taking up to `rescue_candidates` (default 25) references that share **any** k-mer. |
| `all` | Align against every reference in the library. The pre-2026-08-12 behaviour. |
| `none` | Accept the shortlist verdict. |

The default is `kmer` because a read that shares no k-mer with a reference
cannot reach the identity floor anyway, so aligning against it only costs time.

## Why the old path was slow

`assign_sequence` was called once per k-mer size, each time with
`fallback_align_all=True`. Whenever a shortlist failed the identity or coverage
floors, it aligned the read against the entire library — and because `no_match`
did not terminate the outer loop, that whole-library sweep was then repeated for
k=11 and k=9. On the 20k pilot this measured **78 alignments per read**, of which
21,150 of 23,550 came from the brute-force sweep. The restructured path performs
**4.1 alignments per read**.

A second cost was `canonical_unique_kmers`, which normalised and
reverse-complemented every sliding window individually. It now normalises once
per sequence and reverse-complements each unambiguous run once: 4.24 ms → 0.47 ms
for a 1,291 nt read, verified bit-identical over randomised inputs.

## Equivalence evidence

Measured on `runs/pilot_inputs/20260506_lab_biotin_byte_sample_20000.fastq`
(15,126 demultiplexed reads) against
`runs/20260506-lab-biotin-pilot-20k-legacy-optimised`:

- Every assignment status and `reference_ids` value is unchanged, and
  `consensus.fasta`, `consensus.csv.gz`, and `qc.csv.gz` are byte-identical.
- `rescue: kmer` and `rescue: all` were run over all 15,126 reads and agreed on
  every read.
- Recorded evidence differs only on rejected reads. All 740 rows with a changed
  `best_identity` are `no_match` — the *near-miss* reference reported for a read
  that was rejected either way can differ. On assigned rows `best_identity`
  changed in neither direction and `identity_margin` is completely unchanged, so
  the shortlist never missed a reference the exhaustive sweep would have chosen.

`tests/test_assignment_determinism.py` locks in that `rescue="kmer"` and
`rescue="all"` agree on synthetic reads carrying 6% substitutions and indels.

## Determinism

`ReferenceIndex.shortlist` previously accumulated `scores[aliases] += 1.0/n`
while iterating a `frozenset`. Set iteration order for strings depends on
`PYTHONHASHSEED`, which is randomised per process, and floating-point addition is
not associative — so the same read could score differently in different process
workers. This was observed: the first process-backend pilot differed from the
thread-backend pilot on 2 of 15,126 reads.

Scores are now accumulated as exact integer tallies per specificity class and
summed in a canonical order. `tests/test_assignment_determinism.py` runs the
assignment path under three `PYTHONHASHSEED` values and compares scores as
hexadecimal floats; it fails against the pre-fix implementation.

Serial and 8-process runs of the pilot are byte-identical on every artifact.

## Performance and tuning

Assignment is bound by Python-level k-mer work, not by the GIL-releasing edlib
calls, so **threads do not help and processes do**. Steady-state throughput on a
16-logical-CPU Apple Silicon Mac, with per-worker startup removed by solving
`T(n) = startup + n/rate`:

| Backend | Throughput |
| --- | --- |
| serial | 660 reads/s |
| process, 4 | 1,378 reads/s |
| process, 8 | 2,127 reads/s |
| process, 16 | 2,021 reads/s |

Process workers cost about 1 s each to start, because each builds its own
reference indexes rather than receiving them with every task. Short pilots
therefore favour fewer workers than full runs do.

Benchmark on your own machine before changing a profile:

```bash
python scripts/benchmark_assignment.py \
  --config configs/runs/20260506_lab_biotin_pilot20k.yaml \
  --demuxed-reads runs/<run-id>/stages/02_demux/demuxed_reads.jsonl.gz \
  --backend process --jobs 8
```

Add `--rescue all` to compare the exhaustive policy against the default on your
own data before relying on the bounded one.


## Coding QC and the constant regions

Reading-frame and internal-stop QC are driven by constant sequence the user
declares, so they adapt to any construct:

```yaml
qc:
  upstream_constant: ATGCAGCTT     # must begin at the start codon
  downstream_constant: AGTGGATCC...TAA
```

QC assembles `upstream_constant + consensus + downstream_constant`, requires a
whole number of codons, translates frame 0, and fails when a stop appears before
the final codon. Both constants are required together: a frame inferred from one
side is a guess, and stop codons read in a guessed frame are noise rather than
biology. For the same reason `internal_stops` is `not_evaluable` whenever the
frame check fails. Codons containing an ambiguity code translate to `X` rather
than a guessed residue.

**These constants must be kept consistent with the motifs.** `forward_motif`
determines where the extracted insert begins, and `upstream_constant` must cover
exactly the coding sequence between the start codon and that point. On the
20260506 profile the motif runs through the constant `CAGCTT` linker, so
`upstream_constant` is the full `ATGCAGCTT`; before that change the motif stopped
at the ATG and the correct value was `ATG` alone. Changing one without the other
silently shifts the frame.
