# Nanopore2 migration and validation plan

Nanopore2 remains read-only historical evidence. Algorithms are reimplemented from documented
behavior rather than moving its mutable output trees into Nanopore3.

## Implemented foundation

- Strict plain/gzip FASTQ and FASTA readers with stable read UIDs.
- Barcode preflight, best/second evidence, margins and explicit rejection states.
- Mandatory ordered motif extraction.
- Specificity-weighted k-mer shortlist followed by alignment validation.
- Exact duplicate-reference alias sets.
- Deterministic bounded-depth portable consensus with contributor provenance.
- Independent QC criteria with `pass`, `fail` and `not_evaluable`.
- Atomic immutable stages, checksums, safe resume, HTML summary and cross-worker tests.

## Next scientific gates

1. Build Nanopore2-compatible assay profiles from plate/well barcodes, motifs and reference
   manifests; never derive identity from filenames.
2. Simulate ONT substitutions, indels, homopolymers, truncations, orientation changes and known
   two-reference breakpoints over the real reference library.
3. Select thresholds on development data and freeze them before evaluating held-out plates/runs.
4. Benchmark portable consensus against SPOA and MAFFT-based adapters by accuracy, depth and time.
5. Add a split-alignment chimera candidate model. Promote it to a biological call only after its
   false-positive rate against truncations, large deletions and close references is acceptable.
6. Compare clean Nanopore3 runs against archived Nanopore2 outputs by per-read decision—not just
   aggregate assignment yield—and investigate every major discordant category.

## Deferred intentionally

Fragment/FULL_INFO heuristics, automated identity reassignment, advanced overhang/enrichment
statistics, a graphical interface, distributed/cloud executors and the unrelated SUMO analysis
are outside the first validation release. This keeps the initial scientific surface small enough
to test properly.
