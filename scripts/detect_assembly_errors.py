#!/usr/bin/env python3
"""Classify oligo-pool assembly outcomes per read using the Golden Gate design.

For each demultiplexed read this projects the designed fragment boundaries of the
best-matching gene onto the read, then asks which design each fragment actually
came from.  Because two designs can only mis-assemble where they share an
overhang pair, the candidate donors for a slot are exactly the variants carrying
that pair.

Outcomes: ``intact``, ``chimeric`` (fragments from different designs, with the
junction and donors named), ``fragment_missing``, ``fragment_extra``, and
``unresolved`` when the evidence does not separate the possibilities.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import edlib

from nanopore3.assignment import ReferenceIndex, assign_read
from nanopore3.config import load_config
from nanopore3.fragments import FragmentLibrary
from nanopore3.references import read_reference_libraries

_CIGAR = re.compile(r"(\d+)([=XID])")


def project_spans(query: str, gene: str, spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Map gene coordinates onto query coordinates through an NW alignment."""

    result = edlib.align(query, gene, mode="NW", task="path")
    cigar = result.get("cigar") or ""
    # gene position -> query position
    mapping: dict[int, int] = {}
    qi = gi = 0
    for count_text, operation in _CIGAR.findall(cigar):
        count = int(count_text)
        if operation in "=X":
            for _ in range(count):
                mapping[gi] = qi
                qi += 1
                gi += 1
        elif operation == "I":       # present in query, absent from gene
            qi += count
        elif operation == "D":       # present in gene, absent from query
            for _ in range(count):
                mapping[gi] = qi
                gi += 1
    projected = []
    for start, end in spans:
        lo = max(0, min(start, len(gene) - 1))
        hi = max(0, min(end, len(gene)))
        qs = mapping.get(lo, 0)
        qe = mapping.get(hi - 1, len(query) - 1) + 1
        projected.append((qs, max(qs, qe)))
    return projected


def best_variant(
    observed: str, variants: dict[str, tuple[str, ...]]
) -> tuple[tuple[str, ...], float, float]:
    """Return the owners of the closest variant, its identity, and its margin."""

    scored = []
    for sequence, owners in variants.items():
        distance = edlib.align(observed, sequence, mode="NW", task="distance")["editDistance"]
        identity = max(0.0, 1.0 - distance / max(len(observed), len(sequence), 1))
        scored.append((identity, owners))
    scored.sort(key=lambda item: (-item[0], item[1]))
    if not scored:
        return (), 0.0, 0.0
    top_identity, owners = scored[0]
    margin = top_identity - scored[1][0] if len(scored) > 1 else top_identity
    return owners, top_identity, margin


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--demuxed-reads", type=Path, required=True)
    parser.add_argument(
        "--fragments",
        action="append",
        required=True,
        metavar="LIBRARY=CSV",
        help="reference library id and its oPool *_FULL_INFO.csv",
    )
    parser.add_argument("--limit", type=int)
    parser.add_argument("--min-fragment-identity", type=float, default=0.80)
    parser.add_argument("--min-fragment-margin", type=float, default=0.03)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    collection = read_reference_libraries(
        {library: settings.fasta for library, settings in config.reference_sets.items()}
    )
    references = {
        library: {record.id: record.sequence for record in bundle.records}
        for library, bundle in collection.libraries
    }
    indexes = {
        library: tuple(
            ReferenceIndex(
                references[library],
                k=k,
                max_kmer_owners=config.reference_sets[library].max_kmer_owners,
            )
            for k in config.reference_sets[library].kmer_sizes
        )
        for library in collection.ids
    }
    fragment_libraries = {}
    for item in args.fragments:
        library, _, path = item.partition("=")
        fragment_libraries[library] = FragmentLibrary.from_full_info_csv(path)
        print(f"# {library}: {fragment_libraries[library].summary()}")

    outcomes: Counter[str] = Counter()
    donors: Counter[tuple[str, str]] = Counter()
    rows = []
    with gzip.open(args.demuxed_reads, "rt", encoding="utf-8", newline="") as handle:
        for position, line in enumerate(handle):
            if args.limit is not None and position >= args.limit:
                break
            read = json.loads(line)
            library = config.reference_library_id_for_plate(str(read["plate_id"]))
            settings = config.reference_sets[library]
            fragments = fragment_libraries.get(library)
            if fragments is None:
                continue
            call, _ = assign_read(
                read["sequence"], indexes[library],
                left_motif=config.library.forward_motif,
                right_motif=config.library.reverse_motif,
                motif_max_edits=config.library.motif_max_edits,
                top_n=settings.candidate_count,
                min_kmer_score=settings.minimum_kmer_score,
                min_identity=settings.minimum_identity,
                min_query_coverage=settings.minimum_query_coverage,
                min_reference_coverage=settings.minimum_reference_coverage,
                min_identity_margin=settings.minimum_identity_margin,
                rescue=settings.rescue_policy,
            )
            if call.extraction is None or not call.extraction.sequence:
                outcomes["no_insert"] += 1
                continue
            query = call.extraction.sequence
            # Use the best alignment even when it failed the assignment floors:
            # a chimera is expected to fail them, and is still anchored by one parent.
            anchor = call.best.aliases[0] if call.best else None
            if anchor is None:
                outcomes["unresolved"] += 1
                continue
            gene = fragments.gene_for_sequence(references[library][anchor])
            if gene is None:
                outcomes["gene_not_in_design"] += 1
                continue
            if len(gene.fragments) < 2:
                outcomes["single_fragment_gene"] += 1
                continue
            spans = [(f.start, f.end) for f in gene.fragments]
            projected = project_spans(query, gene.sequence, spans)
            slot_owners: list[tuple[str, ...]] = []
            slot_identities: list[float] = []
            resolved = True
            for fragment, (qs, qe) in zip(gene.fragments, projected, strict=True):
                observed = query[qs:qe]
                if len(observed) < 0.5 * len(fragment.sequence):
                    outcomes["fragment_missing"] += 1
                    resolved = False
                    break
                owners, identity, margin = best_variant(
                    observed, fragments.interchangeable(fragment)
                )
                if identity < args.min_fragment_identity or margin < args.min_fragment_margin:
                    resolved = False
                    slot_owners.append(())
                else:
                    slot_owners.append(owners)
                slot_identities.append(identity)
            else:
                if not resolved or any(not o for o in slot_owners):
                    outcomes["unresolved"] += 1
                else:
                    common = set(slot_owners[0])
                    for owners in slot_owners[1:]:
                        common &= set(owners)
                    if common:
                        outcomes["intact"] += 1
                    else:
                        outcomes["chimeric"] += 1
                        pair = (
                            sorted(slot_owners[0])[0],
                            sorted(slot_owners[-1])[0],
                        )
                        donors[pair] += 1
                        rows.append(
                            {
                                "read_uid": read["read_uid"],
                                "plate_id": read["plate_id"],
                                "well_id": read["well_id"],
                                "reference_library_id": library,
                                "assignment_status": call.status,
                                "reported_reference": "|".join(call.reference_ids),
                                "slot_owners": " >> ".join(
                                    sorted(o)[0] for o in slot_owners
                                ),
                                "slot_identities": " ".join(
                                    f"{value:.3f}" for value in slot_identities
                                ),
                            }
                        )
    total = sum(outcomes.values())
    print(f"\n# assembly outcomes over {total} reads")
    for name, count in outcomes.most_common():
        print(f"   {name:24s} {count:6d}  {count/total*100:5.1f}%")
    decided = outcomes["intact"] + outcomes["chimeric"]
    if decided:
        print(f"\n# chimeric among decided reads: {outcomes['chimeric']}/{decided} "
              f"({outcomes['chimeric']/decided*100:.1f}%)")
    if donors:
        print("\n# most frequent donor pairs (5' design >> 3' design)")
        for (left, right), count in donors.most_common(8):
            print(f"   {count:5d}  {left[:34]:36s} >> {right[:34]}")
    if args.out and rows:
        with args.out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n# wrote {len(rows)} chimeric read records to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
