"""Re-run the analysis stages from an existing run's intermediates.

Demultiplexing 2.5 M reads takes most of a run's wall clock and depends only on the
barcodes; changing a consensus threshold and repeating the whole thing wastes it,
and requires the original FASTQ to still be attached.  This starts from a finished
run instead: the stages that would not change are carried over, and everything from
a chosen stage onward is recomputed.

Two properties are kept, because they are what makes a run trustworthy:

* **Runs stay immutable.**  A rerun writes a *new* directory rather than editing the
  source, and records which run it inherited from.  The inherited stages are hard
  linked where the filesystem allows it, so carrying 600 MB of intermediates
  forward costs nothing and still cannot be modified in place.
* **Reuse is decided by the stage fingerprint, not by trust.**  Every stage records
  a fingerprint over its parameters, its input digests and the pipeline version.  A
  stage is inherited only if the *new* configuration produces the same fingerprint
  the source recorded.  Change a barcode setting and the demux fingerprint changes,
  and the rerun refuses to inherit it rather than quietly building on a stage that
  no longer matches the configuration.

That second point is why this is safe without re-reading the original FASTQ: the
inherited stage carries its input digests forward, so provenance still names the
input file and its checksum without the file needing to be present.
"""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

from .errors import Nanopore3Error
from .provenance import (
    StageValidationError,
    require_current_implementation,
    validate_stage_directory,
)

# Every stage, in the order they run. A rerun from a stage recomputes it and all
# that follow.
STAGE_ORDER = (
    "01_ingest",
    "02_demux",
    "03_assignment",
    "03b_chimera",
    "04_consensus",
    "05_qc",
    "06_report",
)

# The stages that read the original FASTQ. Inheriting both of these is what lets a
# rerun proceed with the source data detached.
INPUT_CONSUMING_STAGES = ("01_ingest", "02_demux")


class RerunError(Nanopore3Error, ValueError):
    """The source run cannot serve as a basis for this rerun."""


def stages_before(stage: str) -> tuple[str, ...]:
    """The stages that would be inherited when recomputing from ``stage``."""

    if stage not in STAGE_ORDER:
        raise RerunError(
            f"unknown stage {stage!r}; expected one of {', '.join(STAGE_ORDER)}"
        )
    return STAGE_ORDER[: STAGE_ORDER.index(stage)]


def completed_stages(run_dir: Path) -> tuple[str, ...]:
    """Stages of a run that finished, in pipeline order."""

    stages = run_dir / "stages"
    if not stages.is_dir():
        return ()
    return tuple(
        stage
        for stage in STAGE_ORDER
        if (stages / stage / "_SUCCESS").is_file()
    )


def _link_tree(source: Path, destination: Path) -> None:
    """Copy a stage directory, hard linking files where the filesystem allows.

    A hard link cannot be edited into a different file without breaking the link,
    and stage directories are write-once, so this is safe as well as cheap.  Falls
    back to copying across filesystems.
    """

    destination.mkdir(parents=True, exist_ok=False)
    for item in sorted(source.rglob("*")):
        target = destination / item.relative_to(source)
        if item.is_dir():
            target.mkdir(parents=True, exist_ok=True)
            continue
        try:
            os.link(item, target)
        except OSError:
            shutil.copy2(item, target)


def prepare_rerun(
    source_run: Path,
    destination_run: Path,
    from_stage: str,
    *,
    expected_fingerprints: dict[str, str] | None = None,
    verify_checksums: bool = False,
) -> dict[str, object]:
    """Carry a source run's earlier stages into a new run directory.

    ``expected_fingerprints`` maps stage name to the fingerprint the *new*
    configuration computes for it.  Any inherited stage whose recorded fingerprint
    differs is refused by name: it means the configuration change reaches back into
    a stage this rerun intended to keep, so the intermediates no longer describe
    the configuration being run.
    """

    inherited = stages_before(from_stage)
    if not inherited:
        raise RerunError(
            f"{from_stage} is the first stage, so there is nothing to inherit; "
            "start a normal run instead"
        )
    source_stages = source_run / "stages"
    if not (source_run / "run.json").is_file():
        raise RerunError(f"not a run directory (no run.json): {source_run}")
    source_metadata = json.loads((source_run / "run.json").read_text(encoding="utf-8"))
    require_current_implementation(source_metadata)

    finished = completed_stages(source_run)
    # Some stages are optional - chimera detection only runs when it is enabled -
    # so an absent one is not a broken source run. Inherit what exists and let the
    # pipeline compute the rest; a stage that is absent *and* required will fail
    # on its own inputs, which is a clearer error than anything invented here.
    absent = [stage for stage in inherited if stage not in finished]
    inherited = tuple(stage for stage in inherited if stage in finished)
    if not inherited:
        raise RerunError(
            f"{source_run.name} completed none of the stages before {from_stage}, "
            "so there is nothing to build on"
        )

    manifests: dict[str, dict[str, object]] = {}
    for stage in inherited:
        expected = (expected_fingerprints or {}).get(stage)
        try:
            manifest = validate_stage_directory(
                source_stages / stage,
                expected_fingerprint=expected,
                expected_stage=stage,
                verify_checksums=verify_checksums,
            )
            require_current_implementation(manifest.runtime)
        except StageValidationError as exc:
            if expected is not None:
                raise RerunError(
                    f"cannot inherit {stage}: the current configuration does not "
                    f"produce the stage it recorded ({exc}). Something this rerun "
                    f"meant to keep is affected by the configuration change, so run "
                    f"from {stage} or earlier instead"
                ) from exc
            raise RerunError(f"cannot inherit {stage}: {exc}") from exc
        manifests[stage] = {
            "fingerprint": manifest.fingerprint,
            "input_digests": dict(manifest.input_digests),
        }

    destination_run.mkdir(parents=True, exist_ok=True)
    (destination_run / "stages").mkdir(exist_ok=True)
    for stage in inherited:
        _link_tree(source_stages / stage, destination_run / "stages" / stage)

    return {
        "analysis_implementation": source_metadata["analysis_implementation"],
        "inherited_from": source_metadata.get("run_id", source_run.name),
        "not_inherited": absent,
        "inherited_run_path": str(source_run.resolve()),
        "inherited_stages": list(inherited),
        "recomputed_from": from_stage,
        "stage_fingerprints": manifests,
        # Carried forward so provenance still names the original input and its
        # checksum even though the file was not read again.
        "source_preflight": source_metadata.get("preflight", {}),
    }


__all__ = [
    "INPUT_CONSUMING_STAGES",
    "STAGE_ORDER",
    "RerunError",
    "completed_stages",
    "prepare_rerun",
    "stages_before",
]
