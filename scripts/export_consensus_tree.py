#!/usr/bin/env python3
"""Rebuild the graded per-plate consensus tree for a completed run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from nanopore3.export import write_tree_from_stages


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, help="default: <run>/consensus_by_plate")
    args = parser.parse_args()
    root = args.out or (args.run / "consensus_by_plate")
    summary = write_tree_from_stages(
        args.run / "stages" / "04_consensus",
        args.run / "stages" / "05_qc" / "qc.csv.gz",
        root,
    )
    print(json.dumps(summary["grades"], indent=2))
    print(f"{summary['files_written']} files under {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
