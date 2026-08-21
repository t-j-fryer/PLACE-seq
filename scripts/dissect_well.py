#!/usr/bin/env python3
"""Work out what is actually in one well when its consensus looks wrong.

A `mixed_variants` grade says the reads assigned to one design disagree with each
other beyond the support threshold.  That can mean a mixed culture, a chimeric
clone, or reads forced onto the nearest reference when none of them fits.  This
tells the three apart: it clusters the well's reads without reference to any
design, then asks what each cluster's own consensus matches.

    python scripts/dissect_well.py --run-dir runs/260608-full-length-v5c \
        --plate RP07 --well A12

Clustering is reference-free on purpose.  Asking "which design does this read
look like" is the question that produced the confusing answer; asking "which
reads look like each other" is independent of it.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import edlib  # noqa: E402

from nanopore3 import pipeline  # noqa: E402
from nanopore3.assignment import (
    ReferenceIndex,  # noqa: E402
    extract_insert,  # noqa: E402
)
from nanopore3.chimera import (  # noqa: E402
    classify_origin,
    signature_for,
    synthesise_reference,
)
from nanopore3.clustering import cluster_reads  # noqa: E402
from nanopore3.config import load_config  # noqa: E402
from nanopore3.consensus import ConsensusRead, build_reference_consensus  # noqa: E402
from nanopore3.references import read_reference_libraries  # noqa: E402


def best_match(sequence: str, references: dict[str, str], top: int = 3):
    """The closest references to a sequence, by infix edit distance."""

    scored = []
    for name, reference in references.items():
        result = edlib.align(reference, sequence, mode="HW", task="distance")
        distance = result["editDistance"]
        if distance >= 0:
            scored.append((distance, 1 - distance / len(reference), name, len(reference)))
    scored.sort()
    return scored[:top]


def read_signatures(reads, index, *, window: int, step: int):
    """Per-read positional signature, using the pipeline's own chimera profiler.

    An earlier version of this script scanned windows by hand and aligned each
    whole reference into a 90 nt window, which is meaningless and reported
    junctions that were not there.  The tested profiler does the reverse and is
    the code the chimera stage itself runs.
    """

    return [
        signature_for(read_id, sequence, index, window=window, step=step)
        for read_id, sequence in reads
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, help="default: the run's own config")
    parser.add_argument("--plate", required=True)
    parser.add_argument("--well", required=True)
    parser.add_argument("--min-cluster-size", type=int, default=6)
    parser.add_argument("--window", type=int, default=90)
    parser.add_argument("--step", type=int, default=45)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    run_config = json.loads((args.run_dir / "run.json").read_text(encoding="utf-8"))
    config = load_config(args.config) if args.config else None
    if config is None:
        name = run_config["config"]["run_name"]
        candidates = sorted(Path("configs/runs").glob("*.yaml"))
        for candidate in candidates:
            if load_config(candidate).run_name == name:
                config = load_config(candidate)
                break
    if config is None:
        parser.error("could not find the run's config; pass --config")

    flanks = pipeline.resolve_flanks(config)
    config = pipeline.apply_flanks(config, flanks)
    library = config.reference_library_id_for_plate(args.plate)
    collection = read_reference_libraries(
        {k: s.fasta for k, s in config.reference_sets.items()},
        transforms=pipeline.flank_transforms(flanks),
    )
    full = {r.id: r.sequence for r in collection.get(library).records}
    # Reference-free clustering works on the insert: that is where the designs
    # differ, and 849 constant bases would dominate every distance.
    inserts, motifs = pipeline.insert_view({library: full}, flanks)
    references = inserts[library]
    left, right = motifs.get(library) or (
        config.reference_sets[library].forward_motif or config.library.forward_motif,
        config.reference_sets[library].reverse_motif or config.library.reverse_motif,
    )

    reads: list[tuple[str, str]] = []
    with gzip.open(
        args.run_dir / "stages" / "02_demux" / "demuxed_reads.jsonl.gz", "rt", encoding="utf-8"
    ) as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("plate_id") != args.plate or record.get("well_id") != args.well:
                continue
            extraction = extract_insert(record["sequence"], left, right, max_edits=3)
            if extraction.status == "found" and extraction.sequence:
                reads.append((record["read_uid"], extraction.sequence))

    print(f"{args.plate} {args.well}: {len(reads)} reads with a usable insert region")
    if not reads:
        return 1
    lengths = sorted(len(sequence) for _, sequence in reads)
    print(f"  insert length: min {lengths[0]} median {lengths[len(lengths)//2]} max {lengths[-1]}")

    result = cluster_reads(reads, min_cluster_size=args.min_cluster_size)
    print(f"\nreference-free clustering: {len(result.clusters)} cluster(s), "
          f"{len(result.unassigned_read_ids)} read(s) unassigned"
          + (f" - {result.reason}" if result.reason else ""))
    print(f"  co-varying positions used: {len(result.variant_positions)}")

    report: dict[str, object] = {
        "plate": args.plate,
        "well": args.well,
        "reads": len(reads),
        "clusters": [],
        "unassigned": len(result.unassigned_read_ids),
        "reason": result.reason,
    }
    for cluster in result.clusters:
        hits = best_match(cluster.consensus, references)
        print(f"\n  cluster {cluster.cluster_id}: {cluster.size} reads, "
              f"consensus {len(cluster.consensus)} nt, "
              f"mean mismatches to centre {cluster.mean_mismatches_to_centre:.1f}")
        for distance, identity, name, reference_length in hits:
            print(
                f"      {100 * identity:6.2f}% ({distance:3d} edits "
                f"vs {reference_length} nt)  {name}"
            )
        report["clusters"].append(
            {
                "cluster_id": cluster.cluster_id,
                "reads": cluster.size,
                "consensus_length": len(cluster.consensus),
                "best_matches": [
                    {"reference": n, "identity": round(i, 6), "edit_distance": d}
                    for d, i, n, _ in hits
                ],
            }
        )

    # Reads grouped by what they look like window by window: the question the
    # clustering above cannot answer, because a mosaic read has no single owner.
    index = ReferenceIndex(references, k=15, max_kmer_owners=10)
    signatures = read_signatures(reads, index, window=args.window, step=args.step)
    by_signature: dict[tuple[str, ...], list[str]] = {}
    for signature in signatures:
        by_signature.setdefault(signature.signature, []).append(signature.read_id)
    sequences = dict(reads)
    print(f"\npositional signatures over {len(signatures)} reads: "
          f"{len(by_signature)} distinct")
    groups = []
    for shape, read_ids in sorted(by_signature.items(), key=lambda kv: -len(kv[1])):
        if len(read_ids) < args.min_cluster_size:
            continue
        origin = classify_origin(shape, None) if len(shape) > 1 else "single"
        print(f"\n  {len(read_ids):4d} reads, {len(shape)} segment(s), origin {origin}")
        for name in shape:
            print(f"        {name}")
        entry = {"reads": len(read_ids), "segments": list(shape), "origin": origin}
        if shape:
            scaffold = (
                references.get(shape[0], "")
                if len(shape) == 1
                else synthesise_reference(shape, references, 1, window=args.window, step=args.step)
            )
            members = [ConsensusRead(rid, sequences[rid]) for rid in read_ids[:30]]
            if scaffold and members:
                built = build_reference_consensus(
                    scaffold, members, group_id="dissect", min_depth=1,
                    max_reads=len(members), min_support=0.6, seed=0,
                )
                if built.sequence:
                    hits = best_match(built.sequence, references, top=2)
                    print(f"        -> consensus {len(built.sequence)} nt, "
                          f"{built.ambiguous_bases} ambiguous; closest references:")
                    for distance, identity, name, reference_length in hits:
                        print(f"             {100 * identity:6.2f}% ({distance:3d} edits "
                              f"vs {reference_length} nt)  {name}")
                    entry["consensus_length"] = len(built.sequence)
                    entry["ambiguous_bases"] = built.ambiguous_bases
                    entry["closest"] = [
                        {"reference": n, "identity": round(i, 6), "edit_distance": d}
                        for d, i, n, _ in hits
                    ]
        groups.append(entry)
    report["signature_groups"] = groups

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
