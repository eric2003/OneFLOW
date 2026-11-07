"""Discover OneFLOW projects independently of the built-in case catalogs."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from .script_file_discovery import discover_script_parameter_files


class ProjectType(str, Enum):
    """Known OneFLOW project families inferred from their runtime task."""

    GRID = "grid"
    FLOW = "flow"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class ProjectDefinition:
    """A user-selected project directory and its best-effort classification."""

    name: str
    root: Path
    project_type: ProjectType
    control_file: Path | None
    runtime_file: Path | None
    parameter_files: tuple[Path, ...] = ()
    parameter_file_diagnostics: tuple[str, ...] = ()
    parameter_file_mode: str = "directory"


def _control_file(root: Path) -> Path | None:
    for relative_path in ("script/control.txt", "control.txt"):
        candidate = root / relative_path
        if candidate.is_file():
            return candidate.resolve()
    return None


def _workflow_control_files(root: Path) -> list[Path]:
    """Return case control files in workflow order."""

    workflow_file = root / "workflow.json"
    try:
        workflow = json.loads(workflow_file.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        workflow = {}

    control_files = []
    for step in workflow.get("steps", []):
        for case in step.get("cases", []):
            directory = case.get("directory")
            if not isinstance(directory, str):
                continue
            script = root / directory / "script"
            candidate = script / "control.txt"
            if candidate.is_file() or script.is_dir():
                control_files.append(candidate.resolve())
    if control_files:
        return control_files
    return sorted((path / "control.txt").resolve() for path in root.glob("*/script") if path.is_dir())


def _script_directory(root: Path, control_file: Path | None) -> Path:
    if control_file is not None:
        return control_file.parent
    script = root / "script"
    return script if script.is_dir() else root


def _discover_project_parameters(
    root: Path,
    control_files: list[Path],
    primary_control: Path | None,
) -> tuple[tuple[Path, ...], tuple[str, ...], str]:
    directories = (
        [control.parent for control in control_files]
        if control_files
        else [_script_directory(root, primary_control)]
    )
    discoveries = [
        discover_script_parameter_files(directory) for directory in directories
    ]
    files = tuple(path for discovery in discoveries for path in discovery.files)
    diagnostics = tuple(
        f"{directory}: {message}"
        for directory, discovery in zip(directories, discoveries)
        for message in discovery.diagnostics
    )
    modes = {discovery.mode for discovery in discoveries}
    mode = next(iter(modes)) if len(modes) == 1 else "mixed"
    return files, diagnostics, mode


def _runtime_file(parameter_files: tuple[Path, ...]) -> Path | None:
    ordered = sorted(
        parameter_files,
        key=lambda path: (path.name.casefold() != "simu.txt", path.name),
    )
    for path in ordered:
        try:
            content = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if re.search(r"\bsimutask\s*=", content, flags=re.IGNORECASE):
            return path
    return None

def _detect_type(runtime_file: Path | None) -> ProjectType:
    if runtime_file is None:
        return ProjectType.UNKNOWN

    try:
        content = runtime_file.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ProjectType.UNKNOWN

    task_match = re.search(
        r"""\bsimutask\s*=\s*["']([^"']+)["']""",
        content,
        flags=re.IGNORECASE,
    )
    if task_match is None:
        return ProjectType.UNKNOWN

    task = task_match.group(1).strip().lower()
    if task == "grid":
        return ProjectType.GRID
    if task in {"flow", "cfd", "solver", "solve"}:
        return ProjectType.FLOW
    return ProjectType.UNKNOWN


def open_project(directory: str | Path) -> ProjectDefinition:
    """Open a project and discover its actual parameter files."""

    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Project directory does not exist: {root}")

    is_workflow = (root / "workflow.json").is_file()
    workflow_controls = _workflow_control_files(root) if is_workflow else []
    control_file = _control_file(root)
    if control_file is None and workflow_controls:
        control_file = workflow_controls[0]

    controls_for_parameters = workflow_controls if is_workflow else (
        [control_file] if control_file is not None else []
    )
    parameter_files, diagnostics, mode = _discover_project_parameters(
        root, controls_for_parameters, control_file
    )
    runtime_file = _runtime_file(parameter_files)

    return ProjectDefinition(
        name=root.name or str(root),
        root=root,
        project_type=_detect_type(runtime_file),
        control_file=control_file,
        runtime_file=runtime_file,
        parameter_files=parameter_files,
        parameter_file_diagnostics=diagnostics,
        parameter_file_mode=mode,
    )
