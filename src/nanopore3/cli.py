"""Dependency-light command line interface shared by all supported platforms."""

from __future__ import annotations

import argparse
from importlib.resources import as_file, files
import json
from multiprocessing import freeze_support
from pathlib import Path
import shutil
import sys

from . import __version__
from .config import ConfigError, load_config
from .provenance import StageValidationError
from .pipeline import PipelineError, run_pipeline, validate_inputs
from .runtime import doctor_report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nanopore3",
        description="Conservative, reproducible Nanopore amplicon consensus analysis",
    )
    parser.add_argument("--version", action="version", version=f"nanopore3 {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor", help="report runtime and optional backend capabilities")
    doctor.add_argument("--json", action="store_true", help="emit machine-readable JSON")

    validate = subparsers.add_parser("validate", help="validate configuration, inputs and references")
    validate.add_argument("--config", type=Path, required=True)
    validate.add_argument("--quick", action="store_true", help="do not scan every FASTQ record")

    run = subparsers.add_parser("run", help="execute the portable staged workflow")
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--output", type=Path, help="override the configured run root")
    run.add_argument("--run-id", help="portable directory name for this run")
    run.add_argument("--resume", action="store_true", help="resume only checksum-compatible stages")

    init = subparsers.add_parser("init", help="copy the documented synthetic example")
    init.add_argument("directory", type=Path, nargs="?", default=Path("nanopore3-example"))

    rerun = subparsers.add_parser(
        "rerun",
        help="recompute the analysis stages from a finished run's intermediates",
        description=(
            "Start from a finished run instead of the FASTQ. Stages before --from "
            "are carried over (hard linked, so they cost no disk and stay "
            "immutable) and everything from --from onward is recomputed into a new "
            "run directory. A stage is inherited only if the current configuration "
            "produces the fingerprint that stage recorded, so a configuration "
            "change that reaches back into an inherited stage is refused by name "
            "rather than silently built upon. The original FASTQ need not be "
            "attached when both stages that read it are inherited."
        ),
    )
    rerun.add_argument("--config", type=Path, required=True)
    rerun.add_argument(
        "--from-run", type=Path, required=True, help="the finished run to build on"
    )
    rerun.add_argument(
        "--from",
        dest="from_stage",
        default="04_consensus",
        help="first stage to recompute (default: 04_consensus)",
    )
    rerun.add_argument("--output", type=Path, help="override the configured run root")
    rerun.add_argument("--run-id", help="portable directory name for the new run")
    rerun.add_argument(
        "--verify",
        action="store_true",
        help="re-checksum every inherited artifact instead of trusting its manifest",
    )

    layout = subparsers.add_parser(
        "layout",
        help="write a config's pooling layout out as a CSV, or check one",
    )
    layout.add_argument("--config", type=Path, required=True)
    layout.add_argument(
        "--export",
        type=Path,
        help="write the layout to this CSV so it can be maintained in a spreadsheet "
        "and referenced with compressed_pcr.layout_csv",
    )
    return parser


def _print_doctor(as_json: bool) -> None:
    report = doctor_report()
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print(f"Nanopore3 {__version__}")
    print(f"Python: {report['python']}")
    print(f"Platform: {report['platform']}")
    print(f"CPUs: {report['cpu_count']}")
    print("Portable edlib backend: available")
    for name, details in report["optional_binaries"].items():
        state = details["version"] if details["available"] else "not found (optional)"
        print(f"{name}: {state}")


def _init_example(destination: Path) -> None:
    if destination.exists():
        raise PipelineError(f"destination already exists: {destination}")
    (destination / "configs").mkdir(parents=True)
    (destination / "fixtures").mkdir()
    assets = files("nanopore3").joinpath("example_data")
    for name, target in (
        ("example.yaml", destination / "configs" / "example.yaml"),
        ("references.fasta", destination / "fixtures" / "references.fasta"),
        ("reads.fastq", destination / "fixtures" / "reads.fastq"),
    ):
        with as_file(assets.joinpath(name)) as source:
            shutil.copy2(source, target)
    print(f"Created example at {destination.resolve()}")


def _rerun(
    config,
    *,
    source: Path,
    from_stage: str,
    output_root: Path | None,
    run_id: str | None,
    verify: bool,
) -> Path:
    """Carry a finished run's earlier stages forward and recompute the rest."""

    from datetime import datetime, timezone

    from .provenance import canonical_digest
    from .rerun import RerunError, completed_stages, prepare_rerun, stages_before

    source = source.expanduser().resolve()
    inherited = stages_before(from_stage)
    done = completed_stages(source)
    print(f"{source.name}: completed {', '.join(done) or 'nothing'}")
    print(f"inheriting {', '.join(inherited)}; recomputing from {from_stage}")

    if run_id is None:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_id = f"{stamp}-rerun-{canonical_digest(config.as_dict())[:8]}"
    root = (output_root or config.output_root).expanduser().resolve(strict=False)
    destination = root / run_id
    if destination.exists():
        raise PipelineError(f"run already exists; choose another --run-id: {destination}")

    try:
        record = prepare_rerun(
            source, destination, from_stage, verify_checksums=verify
        )
    except RerunError as exc:
        raise PipelineError(str(exc)) from exc
    print(f"carried {len(record['inherited_stages'])} stage(s) into {destination}")
    if record.get("not_inherited"):
        print(
            f"  ({', '.join(record['not_inherited'])} did not run in the source and "
            "will be computed if the configuration calls for it)"
        )
    try:
        return run_pipeline(
            config,
            output_root=output_root,
            run_id=run_id,
            resume=True,
            inherit=record,
        )
    except StageValidationError as exc:
        # The configuration change reaches back into a stage this rerun meant to
        # keep, so its intermediates no longer describe the configuration being
        # run. Refusing is the point; leaving a half-built run behind is not.
        shutil.rmtree(destination, ignore_errors=True)
        carried = ", ".join(record["inherited_stages"])
        raise PipelineError(
            f"the configuration differs from the one that produced an inherited "
            f"stage ({exc}). Inherited: {carried}. Re-run with --from set to the "
            f"earliest stage your change affects, or run from the FASTQ instead"
        ) from exc
    except BaseException:
        shutil.rmtree(destination, ignore_errors=True)
        raise


def _layout(config_path: Path, export: Path | None) -> None:
    """Show or export the pooling layout, and report it against the references."""

    config = load_config(config_path)
    settings = config.compressed_pcr
    if not settings.enabled:
        print("compressed_pcr is not enabled: every well is its own culture plate.")
        return
    plates = sorted({p for ps in settings.pcr_plates.values() for p in ps})
    blocks = sum(len(by_block) for by_block in settings.blocks.values())
    print(
        f"{len(settings.pcr_plates)} colony-PCR barcode(s), {len(plates)} culture "
        f"plate(s), {blocks} block(s) across {len(settings.blocks)} library/libraries"
    )
    for barcode, sources in sorted(settings.pcr_plates.items()):
        mode = settings.clonality.get(barcode, "unspecified")
        print(f"  {barcode}: {len(sources)} culture plate(s), clonality={mode}")

    from .pipeline import check_pooling_layout
    from .references import read_reference_libraries

    problems: list[str] = []
    try:
        from . import pipeline as _pipeline

        flanks = _pipeline.resolve_flanks(config)
        resolved = _pipeline.apply_flanks(config, flanks)
        collection = read_reference_libraries(
            {k: s.fasta for k, s in resolved.reference_sets.items()},
            transforms=_pipeline.flank_transforms(flanks),
        )
        problems = check_pooling_layout(resolved, collection)
    except (OSError, ValueError) as exc:
        print(f"  (references not checked: {exc})")
    if problems:
        print("\nproblems:")
        for problem in problems:
            print(f"  - {problem}")
    else:
        print("\nlayout agrees with the references.")

    if export:
        from .layout import write_layout_csv

        write_layout_csv(
            export,
            dict(settings.pcr_plates),
            {k: dict(v) for k, v in settings.blocks.items()},
            dict(settings.clonality),
        )
        print(f"\nwrote {export}")
        print("Reference it with:\n  compressed_pcr:\n    enabled: true\n"
              f"    layout_csv: {export.name}")


def main(argv: list[str] | None = None) -> int:
    """Run the CLI and return a process exit code."""

    freeze_support()
    args = _parser().parse_args(argv)
    try:
        if args.command == "doctor":
            _print_doctor(args.json)
        elif args.command == "validate":
            config = load_config(args.config)
            result = validate_inputs(config, scan_fastq=not args.quick)
            print(json.dumps(result, indent=2, sort_keys=True))
        elif args.command == "run":
            config = load_config(args.config)
            path = run_pipeline(
                config,
                output_root=args.output,
                run_id=args.run_id,
                resume=args.resume,
            )
            print(f"Completed run: {path}")
            print(f"Report: {path / 'stages' / '06_report' / 'report.html'}")
        elif args.command == "init":
            _init_example(args.directory)
        elif args.command == "rerun":
            config = load_config(args.config)
            path = _rerun(
                config,
                source=args.from_run,
                from_stage=args.from_stage,
                output_root=args.output,
                run_id=args.run_id,
                verify=args.verify,
            )
            print(f"Completed run: {path}")
            print(f"Report: {path / 'stages' / '06_report' / 'report.html'}")
        elif args.command == "layout":
            _layout(args.config, args.export)
        else:  # pragma: no cover - argparse enforces subcommands
            raise AssertionError(args.command)
        return 0
    except (ConfigError, PipelineError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


__all__ = ["main"]
