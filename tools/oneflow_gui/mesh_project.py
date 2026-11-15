"""User-owned mesh-generation project discovery and lifecycle helpers."""

from __future__ import annotations

import json
import re
import shutil
from pathlib import Path

from .user_data_paths import UserDataPaths

_PROJECT_SCHEMA = "oneflow.mesh-generation-project"
_GENERIC_3D_VALIDATION_CASE = "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra"


def _catalog(repository_root: Path) -> dict:
    value = json.loads(
        (Path(repository_root).resolve() / "grid/case_catalog.json").read_text(
            encoding="utf-8"
        )
    )
    if value.get("schema") != "oneflow.case-catalog":
        raise ValueError("Unsupported OneFLOW Case catalog.")
    return value


def _slug(value: str) -> str:
    result = re.sub(r"[^A-Za-z0-9_-]+", "_", value.strip()).strip("_")
    if not result:
        raise ValueError("Project name must contain letters or numbers.")
    return result[:64]


def mesh_project_templates(repository_root: Path) -> tuple[tuple[str, str], ...]:
    """Return supported closed 3-D STL templates from the Case Catalog."""
    entries = _catalog(repository_root).get("cases", [])
    result = []
    for item in entries:
        if item.get("group") != "STL to Gmsh" or item.get("dimension") != 3:
            continue
        if (item.get("execution") or {}).get("kind") != "python":
            continue
        if item.get("id") not in {
            "stl_to_gmsh_python_sphere_automatic_tetrahedra",
            _GENERIC_3D_VALIDATION_CASE,
        }:
            continue
        result.append((str(item["id"]), str(item.get("title", item["id"]))))
    return tuple(result)


def _entry(repository_root: Path, case_id: str) -> dict:
    entry = next(
        (item for item in _catalog(repository_root)["cases"] if item.get("id") == case_id),
        None,
    )
    if entry is None or entry.get("group") != "STL to Gmsh":
        raise ValueError(f"Unknown STL-to-Gmsh project template: {case_id!r}")
    return entry


def _new_project_directory(name: str) -> Path:
    root = UserDataPaths.discover().mesh_generation
    root.mkdir(parents=True, exist_ok=True)
    base = _slug(name)
    candidate = root / base
    suffix = 2
    while candidate.exists():
        candidate = root / f"{base}_{suffix}"
        suffix += 1
    candidate.mkdir(parents=True)
    return candidate


def _write_manifest(directory: Path, manifest: dict) -> dict:
    path = directory / "project.json"
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return load_mesh_project(directory)


def create_mesh_project_from_template(
    repository_root: Path, template_case_id: str, project_name: str,
) -> dict:
    """Copy a checked-in mesh Case into a writable user-owned project."""
    entry = _entry(repository_root, template_case_id)
    source = Path(repository_root).resolve() / entry["directory"]
    if not (source / "input/model.stl").is_file():
        raise FileNotFoundError(f"Template STL is missing: {source / 'input/model.stl'}")
    if not (source / "script/layout.json").is_file():
        raise FileNotFoundError(f"Template layout is missing: {source / 'script/layout.json'}")
    directory = _new_project_directory(project_name)
    template_directory = directory / "template"
    shutil.copytree(source, template_directory)
    layout = json.loads(
        (template_directory / "script/layout.json").read_text(encoding="utf-8")
    )
    manifest = {
        "schema": _PROJECT_SCHEMA,
        "version": 1,
        "name": directory.name,
        "template_case_id": template_case_id,
        "mesh_case_id": template_case_id,
        "input_stl": "template/input/model.stl",
        "output_mesh": f"results/{directory.name}.msh",
        "mesh_size": float(layout.get("mesh_size", 0.20)),
        "boundary_name": str(layout.get("boundary_name", "wall")),
        "volume_name": str(layout.get("surface_name", "fluid_domain")),
    }
    (directory / "results").mkdir()
    return _write_manifest(directory, manifest)


def create_mesh_project_from_stl(
    repository_root: Path, stl_file: Path, project_name: str | None = None,
) -> dict:
    """Import an STL from the repository database or user disk into a new project."""
    source = Path(stl_file).expanduser().resolve()
    if not source.is_file() or source.suffix.lower() != ".stl":
        raise FileNotFoundError(f"Select an existing STL file: {source}")
    directory = _new_project_directory(project_name or source.stem)
    input_directory = directory / "input"
    input_directory.mkdir()
    shutil.copy2(source, input_directory / "model.stl")
    (directory / "results").mkdir()
    manifest = {
        "schema": _PROJECT_SCHEMA,
        "version": 1,
        "name": directory.name,
        "template_case_id": None,
        "mesh_case_id": _GENERIC_3D_VALIDATION_CASE,
        "input_stl": "input/model.stl",
        "output_mesh": f"results/{directory.name}.msh",
        "mesh_size": 0.20,
        "boundary_name": "wall",
        "volume_name": "fluid_domain",
        "source_stl": str(source),
    }
    return _write_manifest(directory, manifest)


def load_mesh_project(directory: Path) -> dict:
    """Load and validate a prepared project directory."""
    root = Path(directory).expanduser().resolve()
    manifest_path = root / "project.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Not a OneFLOW mesh project (project.json missing): {root}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != _PROJECT_SCHEMA or manifest.get("version") != 1:
        raise ValueError(f"Unsupported mesh project manifest: {manifest_path}")
    stl_relative = Path(str(manifest.get("input_stl", "")))
    if stl_relative.is_absolute() or ".." in stl_relative.parts:
        raise ValueError("Project input_stl must be a project-relative path.")
    stl_path = (root / stl_relative).resolve()
    try:
        stl_path.relative_to(root)
    except ValueError as error:
        raise ValueError("Project input STL must remain inside the project directory.") from error
    if not stl_path.is_file() or stl_path.suffix.lower() != ".stl":
        raise FileNotFoundError(f"Project input STL does not exist: {stl_path}")
    output_value = Path(str(manifest.get("output_mesh", "")))
    output_path = output_value if output_value.is_absolute() else (root / output_value)
    if output_path.suffix.lower() != ".msh":
        raise ValueError("Project output_mesh must use the .msh extension.")
    return {
        **manifest,
        "directory": root,
        "stl_path": stl_path,
        "output_path": output_path.resolve(),
    }


def replace_project_stl(project: dict, stl_file: Path) -> dict:
    """Replace a project's source geometry with a private copy of another STL."""
    root = Path(project["directory"]).resolve()
    source = Path(stl_file).expanduser().resolve()
    if not source.is_file() or source.suffix.lower() != ".stl":
        raise FileNotFoundError(f"Select an existing STL file: {source}")
    manifest_path = root / "project.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    relative = Path(str(manifest.get("input_stl", "")))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Project input_stl must be a project-relative path.")
    destination = (root / relative).resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source != destination:
        shutil.copy2(source, destination)
    manifest["source_stl"] = str(source)
    return _write_manifest(root, manifest)


def save_mesh_project_settings(project: dict, settings: dict) -> dict:
    """Persist editable generation settings without changing the source template."""
    root = Path(project["directory"]).resolve()
    manifest_path = root / "project.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key in ("mesh_size", "mesh_size_min", "mesh_size_max", "boundary_name", "volume_name"):
        if key in settings:
            manifest[key] = settings[key]
    output = Path(settings.get("output_path", project["output_path"])).expanduser().resolve()
    try:
        manifest["output_mesh"] = output.relative_to(root).as_posix()
    except ValueError:
        manifest["output_mesh"] = str(output)
    return _write_manifest(root, manifest)


__all__ = [
    "create_mesh_project_from_stl",
    "create_mesh_project_from_template",
    "load_mesh_project",
    "mesh_project_templates",
    "replace_project_stl",
    "save_mesh_project_settings",
]
