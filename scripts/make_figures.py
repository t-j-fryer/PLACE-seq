#!/usr/bin/env python3
"""Render publication figures for a completed run."""

from __future__ import annotations

import argparse
from pathlib import Path

from nanopore3.config import load_config
from nanopore3.deconvolution import CompressedPcrPlan
from nanopore3.figures import write_all


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, help="a completed run directory")
    parser.add_argument("--config", type=Path, help="run profile, for expected clonality")
    parser.add_argument("--out", type=Path, help="output directory (default: <run>/figures)")
    args = parser.parse_args()

    expected: dict[str, int] = {}
    if args.config is not None:
        config = load_config(args.config)
        if config.compressed_pcr.enabled:
            plan = CompressedPcrPlan(
                config.compressed_pcr.pcr_plates,
                config.compressed_pcr.blocks,
                config.compressed_pcr.clonality,
            )
            for plate in plan.pcr_plates:
                value = plan.expected_clones_per_well(plate)
                if value is not None:
                    expected[plate] = value

    output = args.out or (args.run / "figures")
    for path in write_all(args.run, output, expected):
        print(f"wrote {path} (+ .png)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
