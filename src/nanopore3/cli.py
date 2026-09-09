"""Dependency-light command line interface shared by all supported platforms."""

from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import shutil
import sys
import tempfile
from importlib.resources import as_file, files
from multiprocessing import freeze_support
from pathlib import Path

from . import __version__
from .config import load_config
from .errors import Nanopore3Error
from .pipeline import PipelineError, run_pipeline, validate_inputs
from .provenance import StageValidationError
from .runtime import doctor_report


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nanopore3",
        description="PLACE-seq: Plate and Library Assignment, Consensus and Evaluation",
    )
    parser.add_argument(
        "--version", action="version", version=f"PLACE-seq {__version__} (nanopore3)"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser(
        "doctor", help="report runtime and optional backend capabilities"
    )
    doctor.add_argument("--json", action="store_true", help="emit machine-readable JSON")

    validate = subparsers.add_parser(
        "validate", help="validate configuration, inputs and references"
    )
    validate.add_argument("--config", type=Path, required=True)
    validate.add_argument("--demux-only", action="store_true",
                          help="validate demultiplexing without requiring references")
    validate.add_argument("--quick", action="store_true", help="do not scan every FASTQ record")
    validate.add_argument(
        "--quiet", action="store_true", help="suppress preflight progress"
    )

    run = subparsers.add_parser(
        "run", help="execute the portable staged workflow",
        epilog="YAML defaults: parallel.backend=process, jobs=0 (all allocated CPUs), "
        "threads_per_job=1, group_memory_mb=512 (estimated group memory, not total RSS). "
        "consensus.maximum_reads=0 uses all eligible reads. Resume requires the same "
        "installed source identity; legacy or changed-code runs need a fresh run.",
    )
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--demux-only", action="store_true",
                     help="stop after demultiplexing; export plate/well FASTQ without references")
    run.add_argument("--output", type=Path, help="override the configured run root")
    run.add_argument("--run-id", help="portable directory name for this run")
    run.add_argument("--resume", action="store_true",
                     help="reuse stages with matching source identity, config and checksums")
    run.add_argument(
        "--quiet",
        action="store_true",
        help="suppress per-stage progress (progress is reported by default: a long "
        "run that prints nothing is indistinguishable from one that has hung)",
    )

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
            "attached when both stages that read it are inherited. The installed "
            "source identity must match: legacy runs without an identity or runs "
            "made with different code need a fresh run from original inputs."
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
    rerun.add_argument("--quiet", action="store_true", help="suppress per-stage progress")
    rerun.add_argument(
        "--verify",
        action="store_true",
        help="re-checksum every inherited artifact instead of trusting its manifest",
    )

    subsample = subparsers.add_parser(
        "subsample",
        help="write the first N reads of a FASTQ to a new file",
        description=(
            "For trying a configuration, or checking barcode recovery on a fresh "
            "flowcell, before committing to a full run. Reads are taken in file "
            "order rather than at random, so the result is reproducible and costs "
            "one pass over the head of the file rather than over all of it. "
            "Existing outputs and input aliases are refused. Validated output is "
            "published atomically using a hard link; use local disk in Colab, "
            "then copy to Drive. Unsupported filesystems fail safely."
        ),
    )
    subsample.add_argument("--input", type=Path, required=True)
    subsample.add_argument("--output", type=Path, required=True, help="new file; never overwritten")
    subsample.add_argument("--reads", type=int, default=100_000)

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
    print(f"PLACE-seq {__version__} (nanopore3)")
    print(f"Python: {report['python']}")
    print(f"Platform: {report['platform']}")
    print(f"CPUs: {report['available_cpus']} available ({report['cpu_count']} host)")
    print(
        f"Default group memory estimate budget: {report['group_memory_limit_bytes'] / 2**20:g} MiB "
        "(not a total RSS limit)"
    )
    print(f"Source identity: {report['analysis_implementation']['sha256']}")
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


def _configure_progress(quiet: bool) -> None:
    """Send per-stage progress to stderr unless asked not to.

    Nothing configured logging before, so every LOGGER call in the package went
    nowhere and a run reported only its own completion. On a hosted notebook that
    also risks being disconnected for idleness partway through.
    """

    # Configured on this package's logger rather than the root: raising the root
    # to INFO also turns on every dependency, and matplotlib's font machinery
    # alone buries the progress it was meant to reveal.
    logger = logging.getLogger("nanopore3")
    logger.setLevel(logging.WARNING if quiet else logging.INFO)
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter("%(asctime)s %(message)s", "%H:%M:%S"))
        logger.addHandler(handler)
    logger.propagate = False


def _subsample(source: Path, destination: Path, reads: int) -> int:
    """Validate and atomically publish a FASTQ prefix without replacing a file."""

    if reads < 1:
        raise PipelineError("--reads must be at least 1")
    # lexists also rejects dangling symlinks. Existing outputs include all
    # hard-link/symlink aliases of the input, before either file is opened.
    if os.path.lexists(destination) or source.resolve() == destination.resolve():
        raise PipelineError(f"subsample output already exists or aliases the input: {destination}")
    opener = gzip.open if source.suffix == ".gz" else open
    writer = gzip.open if destination.suffix == ".gz" else open
    destination.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=".subsample-", delete=False
    ) as tmp:
        temporary = Path(tmp.name)
    try:
        # Preserve original bytes, including descriptions, plus lines and CRLF.
        with opener(source, "rb") as handle, writer(temporary, "wb") as out:
            while written < reads:
                block = [handle.readline() for _ in range(4)]
                if not block[0]:
                    break
                if not block[0].startswith(b"@"):
                    raise PipelineError(
                        f"{source} is not a FASTQ: record {written + 1} must begin with '@'"
                    )
                if not all(block):
                    raise PipelineError(
                        f"{source} ends mid-record after {written} complete record(s)"
                    )
                header, sequence, plus, quality = (line.rstrip(b"\r\n") for line in block)
                if (
                    not header[1:].strip()
                    or not plus.startswith(b"+")
                    or len(sequence) != len(quality)
                    or any(c < 33 or c > 126 for c in quality)
                    or any(c not in b"ACGTRYSWKMBDHVNacgtryswkmbdhvn.-" for c in sequence)
                ):
                    raise PipelineError(f"{source} is not a FASTQ: malformed record {written + 1}")
                out.writelines(block)
                written += 1
        # Windows FlushFileBuffers requires a handle opened for writing.
        with temporary.open("r+b") as handle:
            os.fsync(handle.fileno())
        # Unlike replace/rename on POSIX, link is an atomic *no-clobber* create.
        # Unsupported filesystems fail safely, retaining any existing output.
        os.link(temporary, destination)
    except EOFError as exc:
        raise PipelineError(f"{source}: truncated gzip input") from exc
    finally:
        temporary.unlink(missing_ok=True)
    return written


def _report_outputs(path: Path) -> None:
    """Name the three things a person actually opens after a run."""

    print(f"Completed run: {path}")
    if (path / "stages/02_demux/report.html").is_file():
        print("  Demux only: reference assignment, consensus and QC were not run.")
        print(f"  report:  {path / 'stages/02_demux/report.html'}")
        print(f"  FASTQ:   {path / 'stages/02_demux/reads'}")
        print(f"  index:   {path / 'stages/02_demux/reads/index.csv'}")
        return
    print(f"  report:  {path / 'stages' / '06_report' / 'report.html'}")
    tree = path / "consensus_by_plate"
    if tree.is_dir():
        print(f"  clones:  {tree}  (one graded FASTA per clone)")
        print(f"  table:   {tree / 'index.csv'}  (every clone, one row each)")


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

    if config.workflow == "demux_only":
        raise PipelineError("demux-only workflows use `run`, not `rerun`; use --resume to recover")

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
            _configure_progress(args.quiet)
            config = load_config(args.config, demux_only=args.demux_only)
            result = validate_inputs(config, scan_fastq=not args.quick)
            print(json.dumps(result, indent=2, sort_keys=True))
        elif args.command == "run":
            _configure_progress(args.quiet)
            config = load_config(args.config, demux_only=args.demux_only)
            path = run_pipeline(
                config,
                output_root=args.output,
                run_id=args.run_id,
                resume=args.resume,
            )
            _report_outputs(path)
        elif args.command == "init":
            _init_example(args.directory)
        elif args.command == "rerun":
            _configure_progress(args.quiet)
            config = load_config(args.config)
            path = _rerun(
                config,
                source=args.from_run,
                from_stage=args.from_stage,
                output_root=args.output,
                run_id=args.run_id,
                verify=args.verify,
            )
            _report_outputs(path)
        elif args.command == "subsample":
            written = _subsample(args.input, args.output, args.reads)
            print(f"wrote {written:,} read(s) to {args.output}")
        elif args.command == "layout":
            _layout(args.config, args.export)
        else:  # pragma: no cover - argparse enforces subcommands
            raise AssertionError(args.command)
        return 0
    except (Nanopore3Error, OSError) as exc:
        # One base rather than a list of types: a hand-written tuple silently let
        # two error classes through as tracebacks, and would have let through the
        # next one added. Anything that is not a Nanopore3Error is a bug in this
        # package, not a mistake by its user, and should show its traceback.
        print(f"error: {exc}", file=sys.stderr)
        return 2


__all__ = ["main"]
