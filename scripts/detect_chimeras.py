#!/usr/bin/env python3
"""Find chimeric clones in a completed run and build their consensus sequences.

A chimeric molecule is a real clone in a real well.  Assignment cannot describe
one, so its reads are currently lost; this recovers them by matching each read's
windows to the reference library independently, grouping reads that share a
positional signature, and building a consensus per group.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from nanopore3.assignment import ReferenceIndex, extract_insert
from nanopore3.chimera import (
    DEFAULT_STEP,
    DEFAULT_WINDOW,
    group_chimeras,
    signature_for,
    synthesise_reference,
)
from nanopore3.clustering import cluster_reads
from nanopore3.config import load_config
from nanopore3.consensus import ConsensusRead, build_reference_consensus
from nanopore3.deconvolution import block_from_reference_id
from nanopore3.export import normalized_well, safe_name
from nanopore3.references import read_reference_libraries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--plate", action="append", help="limit to these plate barcodes")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    parser.add_argument("--step", type=int, default=DEFAULT_STEP)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    collection = read_reference_libraries(
        {lib: s.fasta for lib, s in config.reference_sets.items()}
    )
    references = {
        lib: {r.id: r.sequence for r in bundle.records}
        for lib, bundle in collection.libraries
    }
    indexes = {
        lib: ReferenceIndex(
            refs, k=15, max_kmer_owners=config.reference_sets[lib].max_kmer_owners
        )
        for lib, refs in references.items()
    }

    # Two designs of the same assembly block reach one colony-PCR well only by
    # being the same colony, so a same-block chimera is a genuine assembled
    # clone. Designs from different blocks arrive as separate colonies pooled
    # into the shared PCR well, so a chimera between them can only have formed
    # in the tube. The distinction is the whole difference between a result and
    # an artefact, so it is recorded per call rather than left to the reader.
    block_pattern = re.compile(config.compressed_pcr.block_pattern)

    def origin(signature) -> str:
        blocks = {block_from_reference_id(p, block_pattern) for p in signature}
        if None in blocks:
            return "unknown"
        return "assembly (same block)" if len(blocks) == 1 else "pcr (different blocks)"

    wanted = set(args.plate) if args.plate else None
    by_well: dict[tuple[str, str, str], list[tuple[str, str]]] = defaultdict(list)
    path = args.run / "stages" / "02_demux" / "demuxed_reads.jsonl.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for line in handle:
            read = json.loads(line)
            plate = str(read["plate_id"])
            if wanted and plate not in wanted:
                continue
            try:
                library = config.reference_library_id_for_plate(plate)
            except KeyError:
                continue
            by_well[(plate, str(read["well_id"]), library)].append(
                (read["read_uid"], read["sequence"])
            )

    root = args.out or (args.run / "chimeras")
    root.mkdir(parents=True, exist_ok=True)
    totals: Counter[str] = Counter()
    rows: list[dict[str, object]] = []

    for (plate, well, library), reads in sorted(by_well.items()):
        settings = config.reference_sets[library]
        index = indexes[library]
        signatures = []
        inserts: dict[str, str] = {}
        for read_id, sequence in reads:
            extraction = extract_insert(
                sequence,
                settings.forward_motif or config.library.forward_motif,
                settings.reverse_motif or config.library.reverse_motif,
                max_edits=(
                    config.library.motif_max_edits
                    if settings.motif_max_edits is None
                    else settings.motif_max_edits
                ),
            )
            if extraction.status != "found" or not extraction.sequence:
                continue
            inserts[read_id] = extraction.sequence
            signatures.append(
                signature_for(read_id, extraction.sequence, index,
                              window=args.window, step=args.step)
            )
        if not signatures:
            continue
        groups, outcome = group_chimeras(
            signatures, minimum_depth=config.consensus.minimum_depth
        )
        totals.update(outcome)
        for group in groups:
            members = [
                ConsensusRead(rid, inserts[rid]) for rid in group.read_ids if rid in inserts
            ]
            if len(members) > config.consensus.maximum_reads:
                members = members[: config.consensus.maximum_reads]
            scaffold = synthesise_reference(
                group.signature, references[library], group.junction_window,
                window=args.window, step=args.step,
            )
            if scaffold:
                result = build_reference_consensus(
                    scaffold, members, group_id=group.label,
                    min_depth=1, max_reads=len(members),
                    min_support=config.consensus.minimum_support,
                )
                sequence, how = result.sequence, "spliced-parent scaffold"
            else:
                # Three or more segments: the junctions are not located well
                # enough to trust a spliced scaffold, so build reference-free.
                clustered = cluster_reads(
                    [(m.read_uid, m.sequence) for m in members],
                    min_cluster_size=1,
                    min_support=config.consensus.minimum_support,
                )
                sequence = clustered.clusters[0].consensus if clustered.clusters else ""
                how = "reference-free"
            if not sequence:
                continue
            directory = root / safe_name(plate) / normalized_well(well)
            directory.mkdir(parents=True, exist_ok=True)
            name = "__".join(safe_name(p, limit=40) for p in group.signature)
            file = directory / f"{safe_name(plate)}_{normalized_well(well)}__{name}__chimera.fasta"
            header = (
                f">{safe_name(plate)}_{normalized_well(well)}_chimera "
                f"parents={'|'.join(group.signature)} reads={group.size} "
                f"origin={origin(group.signature)!r} "
                f"junction_window={group.junction_window} scaffold={how}"
            )
            wrapped = "\n".join(sequence[i:i+60] for i in range(0, len(sequence), 60))
            file.write_text(f"{header}\n{wrapped}\n", encoding="ascii")
            rows.append({
                "plate_id": plate, "well_id": well, "reference_library_id": library,
                "parents": " >> ".join(group.signature), "n_parents": len(group.signature),
                "reads": group.size, "junction_window": group.junction_window,
                "likely_origin": origin(group.signature),
                "scaffold": how, "length": len(sequence),
                "file": str(file.relative_to(root)),
            })
            totals["consensus_written"] += 1
            totals[origin(group.signature)] += 1

    if rows:
        with (root / "chimeras.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (root / "summary.json").write_text(
        json.dumps(dict(sorted(totals.items())), indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(dict(sorted(totals.items())), indent=2))
    print(f"{len(rows)} chimeric consensus sequences written under {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
