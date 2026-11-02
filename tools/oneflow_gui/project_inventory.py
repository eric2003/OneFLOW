"""File inventory and run-to-run change detection for OneFLOW projects."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path


_IGNORED_DIRECTORY_NAMES = {".git", "__pycache__", ".oneflow-runs"}


@dataclass(frozen=True, slots=True)
class ProjectFileEntry:
    """One file or directory visible in a project inventory."""

    relative_path: str
    kind: str
    size_bytes: int
    category: str
    purpose: str
    digest: str = ""


@dataclass(frozen=True, slots=True)
class FileChange:
    """A project entry and its comparison status."""

    entry: ProjectFileEntry
    status: str


def scan_project(root: Path | str) -> dict[str, ProjectFileEntry]:
    """Inventory a project tree without internal VCS, cache, or run-staging files."""
    base = Path(root).expanduser().resolve()
    entries: dict[str, ProjectFileEntry] = {}
    if not base.is_dir():
        return entries

    for path in sorted(base.rglob("*"), key=lambda item: item.as_posix().casefold()):
        relative_parts = path.relative_to(base).parts
        if any(part in _IGNORED_DIRECTORY_NAMES for part in relative_parts):
            continue
        relative = path.relative_to(base).as_posix()
        if path.is_dir():
            entries[relative] = ProjectFileEntry(
                relative, "directory", 0, "Directory", directory_purpose(relative)
            )
            continue
        if not path.is_file():
            continue
        digest = _file_digest(path)
        entries[relative] = ProjectFileEntry(
            relative, "file", path.stat().st_size, classify_category(relative),
            describe_purpose(relative), digest
        )
    return entries


def compare_inventory(
    before: dict[str, ProjectFileEntry],
    after: dict[str, ProjectFileEntry],
) -> tuple[FileChange, ...]:
    """Compare two inventories and return stable path-ordered changes."""
    changes: list[FileChange] = []
    for relative in sorted(set(before) | set(after), key=str.casefold):
        old = before.get(relative)
        new = after.get(relative)
        if old is None:
            changes.append(FileChange(new, "added"))
        elif new is None:
            changes.append(FileChange(old, "deleted"))
        elif old.kind != new.kind or old.size_bytes != new.size_bytes or old.digest != new.digest:
            changes.append(FileChange(new, "modified"))
        else:
            changes.append(FileChange(new, "unchanged"))
    return tuple(changes)


def scan_workflow_run_outputs(
    project_root: Path | str, run_root: Path | str
) -> tuple[ProjectFileEntry, ...]:
    """Return generated or changed files in staged cases, excluding copied inputs."""
    project = Path(project_root).expanduser().resolve()
    run = Path(run_root).expanduser().resolve()
    cases_root = run / "cases"
    if not cases_root.is_dir():
        return ()

    import json

    source_directories: dict[str, Path] = {}
    workflow_file = project / "workflow.json"
    try:
        workflow = json.loads(workflow_file.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        workflow = {}
    for step in workflow.get("steps", []) if isinstance(workflow, dict) else []:
        if not isinstance(step, dict):
            continue
        for case in step.get("cases", []):
            if isinstance(case, dict) and isinstance(case.get("id"), str) and isinstance(case.get("directory"), str):
                source_directories[case["id"]] = (project / case["directory"]).resolve()

    found: list[ProjectFileEntry] = []
    for case_directory in sorted(cases_root.iterdir(), key=lambda item: item.name.casefold()):
        if not case_directory.is_dir():
            continue
        source = source_directories.get(case_directory.name)
        for path in sorted(case_directory.rglob("*"), key=lambda item: item.as_posix().casefold()):
            if not path.is_file():
                continue
            relative_case_path = path.relative_to(case_directory)
            source_file = source / relative_case_path if source is not None else None
            if source_file is not None and source_file.is_file():
                if source_file.stat().st_size == path.stat().st_size and _file_digest(source_file) == _file_digest(path):
                    continue
                status_prefix = "Modified run file"
            else:
                status_prefix = "Generated run file"
            display_path = Path(".oneflow-runs") / run.name / "cases" / case_directory.name / relative_case_path
            relative = display_path.as_posix()
            found.append(ProjectFileEntry(
                relative, "file", path.stat().st_size, classify_category(relative),
                f"{status_prefix}; {describe_purpose(relative)}", _file_digest(path)
            ))
    return tuple(found)


def classify_category(relative_path: str) -> str:
    """Classify a file using stable path and extension conventions."""
    path = Path(relative_path)
    name = path.name.casefold()
    suffix = path.suffix.casefold()
    parts = {part.casefold() for part in path.parts}
    if name == "workflow.json":
        return "Workflow definition"
    if name == "control.txt" or name in {"simu.txt", "dimension.txt", "grid.txt", "cfd.txt"}:
        return "Configuration"
    if "script" in parts or suffix in {".ini", ".cfg", ".toml", ".yaml", ".yml"}:
        return "Configuration"
    if "grid" in parts or suffix in {".grd", ".ofl", ".mesh", ".msh"}:
        return "Mesh data"
    if "results" in parts or suffix in {".plt", ".vtk", ".vtu", ".dat", ".csv"}:
        return "Results / visualization"
    if suffix in {".cpp", ".h", ".hpp", ".c", ".f", ".f90", ".py"}:
        return "Source code"
    if suffix in {".json"}:
        return "Metadata"
    if suffix in {".md", ".rst", ".txt"}:
        return "Documentation / text"
    return "Other"


def describe_purpose(relative_path: str) -> str:
    """Provide a conservative purpose hint based on known OneFLOW conventions."""
    path = Path(relative_path)
    name = path.name.casefold()
    parts = {part.casefold() for part in path.parts}
    known = {
        "workflow.json": "Declares workflow steps, case dependencies, and artifact links.",
        "control.txt": "Selects the runtime parameter files for this case.",
        "simu.txt": "Simulation task and runtime mode settings.",
        "dimension.txt": "Mesh or geometry dimension settings.",
        "grid.txt": "Grid input/output and mesh conversion settings.",
        "cfd.txt": "CFD solver and flow-case settings.",
    }
    if name in known:
        return known[name]
    if path.suffix.casefold() == ".ofl":
        return "OneFLOW mesh file."
    if path.suffix.casefold() in {".grd", ".xyz"}:
        return "Grid coordinates or mesh input."
    if path.suffix.casefold() in {".vtu", ".vtk", ".plt"}:
        return "Visualization or post-processing output."
    if "script" in parts:
        return "Case configuration or runtime input file."
    if "results" in parts:
        return "Generated result or artifact."
    if path.suffix.casefold() in {".cpp", ".h", ".hpp", ".c", ".f", ".f90", ".py"}:
        return "Source code."
    return "Purpose not identified from path or file type."


def directory_purpose(relative_path: str) -> str:
    name = Path(relative_path).name.casefold()
    if name == "script":
        return "Case configuration files."
    if name in {"grid", "mesh"}:
        return "Mesh input or mesh data."
    if name in {"results", "output", "outputs"}:
        return "Generated results and artifacts."
    if name in {"rae2822grid", "rae2822_roe_sa"}:
        return "Workflow case directory."
    return "Project directory."


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
