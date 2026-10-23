"""Discover OneFLOW projects independently of the built-in case catalogs."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path


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


def _control_file(root: Path) -> Path | None:
    for relative_path in ("script/control.txt", "control.txt"):
        candidate = root / relative_path
        if candidate.is_file():
            return candidate.resolve()
    return None


def _parameter_files(control_file: Path | None) -> tuple[Path, ...]:
    if control_file is None:
        return ()

    try:
        content = control_file.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return ()

    files = []
    for value in re.findall(r"""["']([^"']+)["']""", content):
        candidate = (control_file.parent / value).resolve()
        if candidate.is_file():
            files.append(candidate)
    return tuple(files)


def _runtime_file(control_file: Path | None) -> Path | None:
    if control_file is None:
        return None

    try:
        content = control_file.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        return None

    match = re.search(r"""["']([^"']+)["']""", content)
    if match is None:
        return None

    candidate = (control_file.parent / match.group(1)).resolve()
    return candidate if candidate.is_file() else None


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
    if task in {"flow", "cfd", "solver"}:
        return ProjectType.FLOW
    return ProjectType.UNKNOWN


def open_project(directory: str | Path) -> ProjectDefinition:
    """Open any existing directory; catalog registration is not required."""

    root = Path(directory).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Project directory does not exist: {root}")

    control_file = _control_file(root)
    if control_file is None and (root / "workflow.json").is_file():
        case_control_files = sorted(root.glob("*/script/control.txt"))
        if case_control_files:
            control_file = case_control_files[0].resolve()
    runtime_file = _runtime_file(control_file)
    parameter_files = tuple(
        path
        for case_control in (
            [control_file]
            if control_file is not None and not (root / "workflow.json").is_file()
            else sorted(root.glob("*/script/control.txt"))
        )
        for path in _parameter_files(case_control)
    )
    return ProjectDefinition(
        name=root.name or str(root),
        root=root,
        project_type=_detect_type(runtime_file),
        control_file=control_file,
        runtime_file=runtime_file,
        parameter_files=parameter_files,
    )
