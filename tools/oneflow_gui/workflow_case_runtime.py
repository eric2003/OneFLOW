"""Prepare isolated case directories for one workflow graph run."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Mapping

from .workflow_graph_runner import WorkflowGraphRunResult
from .workflow_plan import WorkflowPlanStep
from .workflow_runtime_overrides import apply_runtime_overrides


def prepare_case_directory(
    step: WorkflowPlanStep,
    run_root: Path | str,
    runtime_directories: Mapping[str, Path | str],
) -> Path:
    """Copy a case into the run area and apply its declared artifact inputs.

    runtime_directories maps producer case IDs to their staged run directories.
    The caller should invoke this only after dependencies have completed.
    """
    root = Path(run_root).expanduser().resolve()
    target = root / "cases" / step.case_id
    source = step.directory.resolve()

    if target == source or source in target.parents or target in source.parents:
        raise ValueError("Run directory must be separate from the source case tree.")
    if not source.is_dir():
        raise FileNotFoundError(f"Case directory not found: {source}")
    if target.exists():
        raise FileExistsError(f"Staged case directory already exists: {target}")

    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(source, target)

    overrides_by_file: dict[Path, dict[str, str]] = {}
    for binding in step.inputs:
        producer_directory = runtime_directories.get(binding.producer_id)
        if producer_directory is None:
            raise ValueError(
                f"Missing runtime directory for producer '{binding.producer_id}'."
            )
        parameter_file = target / binding.parameter_file
        overrides = overrides_by_file.setdefault(parameter_file, {})
        overrides["sourceGridCaseDir"] = str(Path(producer_directory).resolve())
        overrides[binding.consumer_setting] = binding.relative_path.as_posix()

    for parameter_file, overrides in overrides_by_file.items():
        try:
            parameter_file.resolve().relative_to(target)
        except ValueError as error:
            raise ValueError("Input parameter file must stay within the staged case.") from error
        apply_runtime_overrides(parameter_file, overrides)

    return target
