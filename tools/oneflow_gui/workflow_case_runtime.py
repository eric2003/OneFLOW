"""Prepare isolated case directories for one workflow graph run."""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Mapping

from .path_utils import is_path_within, normalize_path, paths_equal
from .workflow_plan import WorkflowPlanStep
from .workflow_runtime_overrides import apply_runtime_overrides


def prepare_case_directory(
    step: WorkflowPlanStep,
    run_root: Path | str,
    runtime_directories: Mapping[str, Path | str],
) -> Path:
    """Copy a case and apply artifact bindings to the staged parameter files."""
    root = normalize_path(run_root)
    source = normalize_path(step.directory)
    runtime_directory = step.runtime_directory or Path("cases") / step.case_id
    target = normalize_path(root / runtime_directory)
    if not is_path_within(target, root):
        raise ValueError("Staged case directory must stay within the run workspace.")
    if (
        paths_equal(source, target)
        or is_path_within(source, target)
        or is_path_within(target, source)
    ):
        raise ValueError("Run directory must be separate from the source case tree.")
    if not source.is_dir():
        raise FileNotFoundError(f"Case directory not found: {source}")
    if target.exists():
        if not target.is_dir():
            raise FileExistsError(f"Staged case path is not a directory: {target}")
        # Reuse the stable validation workspace. Overlay the case template so
        # existing output directories (which may be open in a Windows viewer)
        # are not recursively deleted.
        shutil.copytree(source, target, dirs_exist_ok=True)
    else:
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, target)

    overrides_by_file: dict[Path, dict[str, str]] = {}
    for output in step.outputs:
        output_path = target / output.relative_path
        if output_path.is_file():
            output_path.unlink()
        elif output_path.exists() and not output_path.is_dir():
            raise ValueError(f"Declared output path has an unsupported type: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output.parameter_file is not None and output.producer_setting is not None:
            parameter_file = target / output.parameter_file
            overrides = overrides_by_file.setdefault(parameter_file, {})
            overrides[output.producer_setting] = output.relative_path.as_posix()

    for binding in step.inputs:
        producer = runtime_directories.get(binding.producer_id)
        if producer is None:
            raise ValueError(f"Missing runtime directory for producer '{binding.producer_id}'.")
        parameter_file = target / binding.parameter_file
        overrides = overrides_by_file.setdefault(parameter_file, {})
        overrides["sourceGridCaseDir"] = str(Path(producer).resolve())
        overrides[binding.consumer_setting] = binding.relative_path.as_posix()

    for parameter_file, overrides in overrides_by_file.items():
        try:
            parameter_file.resolve().relative_to(target)
        except ValueError as error:
            raise ValueError("Runtime parameter file must stay within the staged case.") from error
        apply_runtime_overrides(parameter_file, overrides)
    return target


def prepare_case_directory_in_place(
    step: WorkflowPlanStep,
    runtime_directories: Mapping[str, Path | str],
) -> Path:
    """Prepare a case in its project directory instead of creating a staged copy."""
    target = step.directory.resolve()
    if not target.is_dir():
        raise FileNotFoundError(f"Case directory not found: {target}")

    overrides_by_file: dict[Path, dict[str, str]] = {}
    for output in step.outputs:
        output_path = (target / output.relative_path).resolve()
        try:
            output_path.relative_to(target)
        except ValueError as error:
            raise ValueError("Declared output path must stay within the case.") from error
        if output_path.is_file():
            output_path.unlink()
        elif output_path.exists():
            raise ValueError(f"Declared output path is not a file: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output.parameter_file is not None and output.producer_setting is not None:
            parameter_file = (target / output.parameter_file).resolve()
            overrides = overrides_by_file.setdefault(parameter_file, {})
            overrides[output.producer_setting] = output.relative_path.as_posix()

    for binding in step.inputs:
        producer = runtime_directories.get(binding.producer_id)
        if producer is None:
            raise ValueError(f"Missing runtime directory for producer '{binding.producer_id}'.")
        parameter_file = (target / binding.parameter_file).resolve()
        overrides = overrides_by_file.setdefault(parameter_file, {})
        overrides["sourceGridCaseDir"] = str(Path(producer).resolve())
        overrides[binding.consumer_setting] = binding.relative_path.as_posix()

    for parameter_file, overrides in overrides_by_file.items():
        try:
            parameter_file.relative_to(target)
        except ValueError as error:
            raise ValueError("Runtime parameter file must stay within the case.") from error
        apply_runtime_overrides(parameter_file, overrides)
    return target
