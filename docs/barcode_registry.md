# Barcode registry and assay routing

Nanopore3 keeps barcode identity separate from experimental meaning. A reverse
primer's full sequence identifies a primer family and vector context; a run
profile separately maps a called plate barcode to the reference library that was
actually loaded on that plate.

## Canonical files

- `configs/barcodes/reverse_primer_families.csv` contains 28 full reverse
  primers: `RP01`–`RP12` (pET amplicon), `RWV1`–`RWV8` (whole pET vector), and
  `RBNPPh1`–`RBNPPh8` (phagemid).
- `configs/barcodes/well_forward_primers_a1_h12.csv` contains the exact 96 full
  well primers recovered from `WDM_INLINE_BARCODES` in the Nanopore2 notebook.

Run YAML selects one CSV `family_id`; families are never pooled automatically.
CSV contents and SHA-256 are expanded into the resolved run configuration.

## Known legacy conflicts

- Well tags map column-major: A1=NB01, B1=NB02, …, H1=NB08, A2=NB09, …,
  H12=NB96.
- The exact legacy B1 full primer is 59 nt and lacks the first base expected from
  NB02; it is preserved verbatim pending evidence from reads.
- Old `RPPh01`–`RPPh08` names do not numerically match newer
  `RBNPPh1`–`RBNPPh8` names. They must be related by exact sequence, never by
  number.
- Short `RB01`–`RB07` entries in LevSeq plate FASTAs are bare internal tags, not
  interchangeable with full RP/RWV/phagemid primers.
- `bo_barcodes_fwd.fasta` is not used as canonical well-primer truth because its
  orientation and numbering shift relative to the 96-well NB set.

## 20260506 LAB BIOTIN profile

The versioned run profile maps RP01–RP04 to AAseq Biotin, RP05 to dTF141/dTF142,
and RP06–RP07 to SUMO LAB. The first six mappings were inferred from the original
demultiplexed tree and consensus filenames; the RP07 mapping was confirmed by the
experiment owner on 2026-08-12. RP08+ remain explicitly unmapped until confirmed.

The optimisation notebook recovered two distinct legacy decision rules. Plate
demultiplexing assigns only when exactly one RP identity is under its threshold;
the active well demultiplexer assigns the unique minimum-distance well and calls
only a best-distance tie ambiguous. The profile makes those policies explicit
and uses the recovered conservative RP setting (`trim=12`, `window=400`,
`max_edits=6`) and rescored well setting (`12/400/7`). See
`demultiplex_optimisation.md` for validation limitations and alternatives.
