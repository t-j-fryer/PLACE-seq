# References: four ways to describe what you expect to find

Nanopore3 compares each consensus to a reference. The input is always an amplicon;
what you supply here says **where in it the designed region sits**, so that errors
in your design are never averaged together with errors in the constant sequence
around it.

That constant sequence may be a couple of hundred bases of expression cassette or
several kilobases of vector. It makes no difference to any of this: "cassette" and
"whole vector" are the same case at different lengths, not two modes.

| you have | use | you get |
|---|---|---|
| the designed inserts, and you only want them scored | **plain** | one accuracy figure per clone |
| the designed inserts, and you know the constant regions | **named flanks** | whole-amplicon consensus, per-region accuracy |
| the designed inserts, plus one assembled example | **template** | the same, without transcribing the constant regions |
| the library already assembled — cassettes or whole vectors | **derive** | the same, with nothing to declare |

## Plain: the designed inserts alone

```yaml
reference_libraries:
  my_library:
    fasta: designs.fasta
```

The consensus spans whatever the primers bound, and the reference is the insert, so
that is what accuracy is measured over. The right choice when you do not care about
the constant sequence — not only when the amplicon happens to be short.

## Named flanks

If the amplicon extends past the insert into constant vector sequence, give that
sequence and the references become whole amplicons:

```yaml
reference_libraries:
  my_library:
    fasta: designs.fasta
    flanks:
      upstream: CATAATCCGCACGCATCTGG...
      downstream: AGTGGATCCGAAACACCAGG...
      anchor_length: 20
```

`upstream` and `downstream` are the complete constant regions **including the
primer binding sites at their outer ends**. The outer `anchor_length` bases of
each become the primer anchors: they bound the reconstructed region and are
excluded from it, so no separate `forward_motif`/`reverse_motif` has to be
restated.

## Template: one example construct

Transcribing two long constant regions correctly is the failure mode this avoids.
Give one construct that already contains one of your inserts, and the constant
regions are found by locating that insert in it:

```yaml
    flanks:
      template: example_construct.fasta   # exactly one record
```

The template may be in either orientation. If no insert aligns to it above 0.90
identity the run stops and says so, rather than deriving flanks from a construct
that is not from this library.

## Derive: the references are already assembled

If your FASTA already holds complete sequences — assembled expression cassettes, or
whole vectors, the usual case when a cloning tool emitted them — declare nothing but
`derive`. **Both are the same case**: the shared backbone is found by comparison, so
its length is irrelevant.

```yaml
reference_libraries:
  my_library:
    fasta: assembled_vectors.fasta
    flanks:
      derive: true
      anchor_length: 20
```

The shared backbone is the longest common prefix and suffix across the set, and
what remains is the insert. No alignment is involved: constructs from one backbone
share exact sequence at both ends.

Checked when it runs:

- **at least two references**, since one cannot hold anything constant
- **at least 60 nt shared at each end**, or it reports that they do not share a
  backbone rather than splitting them somewhere arbitrary
- **a variable region exists**, or it reports that they are identical

Verified on the 260608 library: from 684 assembled constructs it recovers the
135 nt and 754 nt constant regions and both primer anchors byte-identically to the
declared ones.

## Why this matters beyond tidiness

Three things become available only once the insert's position is known:

1. **Per-region accuracy.** An error in a designed insert and an error in the
   vector mean different things, and a single whole-amplicon identity averages
   them. `insert_identity`, `flank_5p_identity` and `flank_3p_identity` are
   reported separately.
2. **Insert-scoped acceptance floors** (`insert_thresholds`, on by default).
   Constant flanks were ~70% of every reference in the 260608 run, so a read
   carrying 250 nt of *foreign* insert still scored 0.93 identity over the whole
   amplicon and cleared floors that had been set when the reference was the insert
   alone. Over the insert it covers 0.53 and fails. Thresholds are only meaningful
   relative to the sequence they are measured over.
3. **Insert-scoped chimera detection.** Positional profiling only discriminates
   where references differ. Run over 1,200 bases of which 850 are identical it
   produces weaker signatures and fewer calls — without failing or warning.

A plain full-length library — assembled constructs with no `derive`, no template
and no named flanks — gets none of the three, and nothing says so. If your
references are whole constructs, use `derive: true`.

## Choosing `anchor_length`

The anchors must be unique within the constant regions, or extraction could bound
the wrong copy. If a run stops with

```
the 3' anchor 'GCCATTGACC' occurs 2 times in the constant regions
```

raise `anchor_length` until it is unique. 20 is the default and is usually enough.
