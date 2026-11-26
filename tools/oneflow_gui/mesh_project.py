"""User-owned mesh-generation project discovery and lifecycle helpers."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path

from .path_utils import (
    is_path_within,
    normalize_path,
    paths_equal,
    relative_path_within,
)
from .user_data_paths import UserDataPaths

_PROJECT_SCHEMA = "oneflow.mesh-generation-project"
_GENERIC_3D_VALIDATION_CASE = "stl_to_gmsh_generic_3d_mesh_handoff"
_LEGACY_GENERIC_3D_VALIDATION_CASE = "stl_to_gmsh_python_prolate_spheroid_automatic_tetrahedra"


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
        if item.get("mesh_project_template") is not True:
            continue
        result.append((str(item["id"]), str(item.get("title", item["id"]))))
    return tuple(result)


def _entry(repository_root: Path, case_id: str) -> dict:
    entry = next(
        (item for item in _catalog(repository_root)["cases"] if item.get("id") == case_id),
        None,
    )
    if (
        entry is None
        or entry.get("group") != "STL to Gmsh"
        or entry.get("mesh_project_template") is not True
    ):
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


def _validate_region_metadata_compatibility(stl_file: Path) -> dict | None:
    """Read optional generated-geometry metadata without discarding region identity."""
    manifest_path = Path(stl_file).with_suffix(".geometry.json")
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(
            f"Cannot read STL geometry manifest: {manifest_path}: {error}"
        ) from error
    if manifest.get("schema") != "oneflow.stl-geometry-manifest":
        raise ValueError(f"Unsupported STL geometry manifest: {manifest_path}")
    return manifest


def _copy_region_manifest(source_stl: Path, destination_stl: Path) -> dict | None:
    """Copy checksum-verified named boundary STL regions beside a project STL."""
    manifest_path = source_stl.parent / "region_manifest.json"
    if not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"Cannot read STL region manifest: {manifest_path}: {error}") from error
    if (
        manifest.get("schema") != "oneflow.stl-surface-regions"
        or manifest.get("version") != 1
        or not isinstance(manifest.get("regions"), dict)
        or set(manifest["regions"]) != {"wall", "inflow", "outflow"}
    ):
        raise ValueError(f"Unsupported STL region manifest: {manifest_path}")
    import hashlib
    if hashlib.sha256(source_stl.read_bytes()).hexdigest() != manifest.get("closed_surface_sha256"):
        raise ValueError("Closed STL checksum does not match region_manifest.json")
    for name in ("wall", "inflow", "outflow"):
        record = manifest["regions"][name]
        relative = Path(str(record.get("file", "")))
        if relative.is_absolute() or len(relative.parts) != 1 or ".." in relative.parts:
            raise ValueError(f"{name} region filename must be a simple relative filename")
        region_source = source_stl.parent / relative
        if not region_source.is_file():
            raise FileNotFoundError(f"Missing {name} region STL: {region_source}")
        if hashlib.sha256(region_source.read_bytes()).hexdigest() != record.get("sha256"):
            raise ValueError(f"{name} region STL checksum does not match its manifest")
        shutil.copy2(region_source, destination_stl.parent / relative.name)
    destination_manifest = destination_stl.parent / "region_manifest.json"
    shutil.copy2(manifest_path, destination_manifest)
    return manifest


def create_mesh_project_from_stl(
    repository_root: Path, stl_file: Path, project_name: str | None = None,
) -> dict:
    """Import an STL from the repository database or user disk into a new project."""
    source = Path(stl_file).expanduser().resolve()
    if not source.is_file() or source.suffix.lower() != ".stl":
        raise FileNotFoundError(f"Select an existing STL file: {source}")
    geometry_manifest = _validate_region_metadata_compatibility(source)
    directory = _new_project_directory(project_name or source.stem)
    input_directory = directory / "input"
    input_directory.mkdir()
    project_stl = input_directory / "model.stl"
    shutil.copy2(source, project_stl)
    if geometry_manifest is not None:
        shutil.copy2(
            source.with_suffix(".geometry.json"),
            project_stl.with_suffix(".geometry.json"),
        )
    region_manifest = _copy_region_manifest(source, project_stl)
    (directory / "results").mkdir()
    manifest = {
        "schema": _PROJECT_SCHEMA,
        "version": 1,
        "name": directory.name,
        "template_case_id": None,
        "mesh_case_id": (
            "stl_to_gmsh_pipe_bend_90_degree"
            if geometry_manifest
            and geometry_manifest.get("geometry") == "pipe_bend_90_degree"
            else _GENERIC_3D_VALIDATION_CASE
        ),
        "input_stl": "input/model.stl",
        "output_mesh": f"results/{directory.name}.msh",
        "mesh_size": 0.20,
        "boundary_name": "wall",
        "volume_name": "fluid_domain",
        "source_stl": str(source),
        "geometry": geometry_manifest.get("geometry") if geometry_manifest else None,
        "boundary_regions": sorted(region_manifest["regions"]) if region_manifest else None,
    }
    return _write_manifest(directory, manifest)


def import_mesh_project_directory(repository_root: Path, source_directory: Path) -> dict:
    """Copy a prepared project directory into user data and adopt its STL input."""
    source = Path(source_directory).expanduser().resolve()
    if not source.is_dir():
        raise NotADirectoryError(f"Select an existing project directory: {source}")
    preferred = source / "input/model.stl"
    if preferred.is_file():
        selected_stl = preferred
    else:
        candidates = sorted(path for path in source.rglob("*.stl") if path.is_file())
        if len(candidates) != 1:
            raise ValueError(
                "Prepared project must contain input/model.stl or exactly one STL file."
            )
        selected_stl = candidates[0]
    directory = _new_project_directory(source.name)
    copied = directory / "source_project"
    shutil.copytree(
        source,
        copied,
        ignore=shutil.ignore_patterns(".git", "build", "results", "__pycache__"),
    )
    selected_stl = normalize_path(selected_stl)
    _validate_region_metadata_compatibility(selected_stl)
    if not is_path_within(selected_stl, source):
        raise ValueError("Prepared project STL must remain inside the selected directory.")
    relative_stl = relative_path_within(selected_stl, source)
    copied_stl = copied / relative_stl
    layout_path = copied / "script/layout.json"
    layout = {}
    if layout_path.is_file():
        layout = json.loads(layout_path.read_text(encoding="utf-8"))
    manifest = {
        "schema": _PROJECT_SCHEMA,
        "version": 1,
        "name": directory.name,
        "template_case_id": None,
        "mesh_case_id": _GENERIC_3D_VALIDATION_CASE,
        "input_stl": (Path("source_project") / relative_stl).as_posix(),
        "output_mesh": f"results/{directory.name}.msh",
        "mesh_size": float(layout.get("mesh_size", 0.20)),
        "boundary_name": str(layout.get("boundary_name", "wall")),
        "volume_name": str(layout.get("surface_name", "fluid_domain")),
        "source_directory": str(source),
    }
    (directory / "results").mkdir()
    return _write_manifest(directory, manifest)



def _project_geometry_lineage(
    root: Path, stl_path: Path, output_path: Path,
) -> dict:
    """Discover project geometry lineage using only checksum-verified records."""
    source_digest = hashlib.sha256(stl_path.read_bytes()).hexdigest()
    try:
        source_label = stl_path.relative_to(root).as_posix()
    except ValueError:
        source_label = str(stl_path)
    lineage = {
        "provenance": {
            "kind": "unknown",
            "notes": "No checksum-verified geometry asset manifest is available.",
        },
        "source_representation": {
            "format": "STL", "path": source_label, "sha256": source_digest,
        },
    }
    manifest_path = stl_path.with_suffix(".geometry.json")
    if manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if (
                isinstance(manifest, dict)
                and manifest.get("schema") == "oneflow.stl-geometry-manifest"
                and manifest.get("stl_sha256") == source_digest
            ):
                lineage.update({
                    "asset_id": manifest.get("asset_id"),
                    "definition_sha256": manifest.get("definition_sha256"),
                    "generator": manifest.get("generator", {}),
                    "provenance": {
                        "kind": "generated",
                        "notes": "The project STL matches its geometry manifest checksum.",
                    },
                })
                patches = manifest.get("patch_triangle_indices")
                surface = manifest.get("surface") or {}
                triangle_count = surface.get("triangle_count")
                if (
                    isinstance(patches, dict)
                    and patches
                    and all(
                        isinstance(name, str)
                        and name
                        and isinstance(indices, list)
                        and all(
                            isinstance(index, int)
                            and not isinstance(index, bool)
                            and index >= 0
                            and (
                                not isinstance(triangle_count, int)
                                or index < triangle_count
                            )
                            for index in indices
                        )
                        for name, indices in patches.items()
                    )
                ):
                    lineage["boundary_patches"] = sorted(patches)
            else:
                lineage["manifest_status"] = "invalid_or_checksum_mismatch"
        except (OSError, json.JSONDecodeError):
            lineage["manifest_status"] = "unreadable"
    if output_path.is_file():
        lineage["derived_mesh"] = {
            "format": "Gmsh MSH",
            "path": output_path.relative_to(root).as_posix()
                if is_path_within(output_path, root) else str(output_path),
            "sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
            "parent": "source_representation",
        }

    downstream = []
    validation_root = root / "validation"
    if validation_root.is_dir():
        suffix = ".geometry-lineage.json"
        for sidecar in sorted(validation_root.rglob("*" + suffix)):
            artifact = sidecar.with_name(sidecar.name[:-len(suffix)])
            if not artifact.is_file():
                continue
            try:
                record = json.loads(sidecar.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            metadata = record.get("artifact") if isinstance(record, dict) else None
            geometry_lineage = record.get("geometry_lineage") if isinstance(record, dict) else None
            if (
                not isinstance(record, dict)
                or record.get("schema") != "oneflow.geometry-artifact-lineage"
                or not isinstance(metadata, dict)
                or metadata.get("file_name") != artifact.name
                or metadata.get("sha256") != hashlib.sha256(artifact.read_bytes()).hexdigest()
                or not isinstance(geometry_lineage, dict)
            ):
                continue
            source_record = geometry_lineage.get("source_representation") or {}
            if source_record and source_record.get("sha256") != source_digest:
                continue
            mesh_record = geometry_lineage.get("derived_mesh") or {}
            if output_path.is_file() and mesh_record:
                if mesh_record.get("sha256") != hashlib.sha256(output_path.read_bytes()).hexdigest():
                    continue
            parent = record.get("parent_artifact")
            if isinstance(parent, dict) and parent.get("path"):
                parent_relative = Path(str(parent["path"]))
                if parent_relative.is_absolute() or ".." in parent_relative.parts:
                    continue
                parent_path = normalize_path(validation_root / parent_relative)
                if (
                    not is_path_within(parent_path, validation_root)
                    or not parent_path.is_file()
                    or parent.get("sha256") != hashlib.sha256(parent_path.read_bytes()).hexdigest()
                ):
                    continue
            downstream.append({
                "role": str(metadata.get("role", "Derived artifact")),
                "format": str(metadata.get("format", artifact.suffix.lstrip(".").upper())),
                "path": artifact.relative_to(root).as_posix()
                    if is_path_within(artifact, root) else str(artifact),
                "sha256": metadata["sha256"],
                "parent": record.get("parent_artifact"),
                "sidecar": sidecar.relative_to(root).as_posix(),
            })
    if downstream:
        lineage["downstream_artifacts"] = downstream
    return lineage


def load_mesh_project(directory: Path) -> dict:
    """Load and validate a prepared project directory."""
    root = Path(directory).expanduser().resolve()
    manifest_path = root / "project.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Not a OneFLOW mesh project (project.json missing): {root}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != _PROJECT_SCHEMA or manifest.get("version") != 1:
        raise ValueError(f"Unsupported mesh project manifest: {manifest_path}")
    if (
        manifest.get("template_case_id") is None
        and manifest.get("mesh_case_id") == _LEGACY_GENERIC_3D_VALIDATION_CASE
    ):
        # Earlier imported-STL projects borrowed the prolate-spheroid CFD Case
        # as a generic fallback. Migrate them to the geometry-neutral handoff
        # route so arbitrary meshes are not sent to an incompatible CFD model.
        manifest["mesh_case_id"] = _GENERIC_3D_VALIDATION_CASE
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
            newline="\n",
        )
    stl_relative = Path(str(manifest.get("input_stl", "")))
    if stl_relative.is_absolute() or ".." in stl_relative.parts:
        raise ValueError("Project input_stl must be a project-relative path.")
    stl_path = normalize_path(root / stl_relative)
    if not is_path_within(stl_path, root):
        raise ValueError("Project input STL must remain inside the project directory.")
    if not stl_path.is_file() or stl_path.suffix.lower() != ".stl":
        raise FileNotFoundError(f"Project input STL does not exist: {stl_path}")
    output_value = Path(str(manifest.get("output_mesh", "")))
    if output_value.is_absolute():
        output_path = output_value
    else:
        if ".." in output_value.parts:
            raise ValueError("Project output_mesh must not escape the project directory.")
        output_path = normalize_path(root / output_value)
        if not is_path_within(output_path, root):
            raise ValueError("Project output_mesh must remain inside the project directory.")
    if output_path.suffix.lower() != ".msh":
        raise ValueError("Project output_mesh must use the .msh extension.")
    return {
        **manifest,
        "directory": root,
        "stl_path": stl_path,
        "output_path": output_path.resolve(),
        "geometry_lineage": _project_geometry_lineage(root, stl_path, output_path.resolve()),
    }


def replace_project_stl(project: dict, stl_file: Path) -> dict:
    """Replace a project's source geometry with a private copy of another STL."""
    root = Path(project["directory"]).resolve()
    source = Path(stl_file).expanduser().resolve()
    if not source.is_file() or source.suffix.lower() != ".stl":
        raise FileNotFoundError(f"Select an existing STL file: {source}")
    geometry_manifest = _validate_region_metadata_compatibility(source)
    manifest_path = root / "project.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    relative = Path(str(manifest.get("input_stl", "")))
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Project input_stl must be a project-relative path.")
    destination = normalize_path(root / relative)
    if not is_path_within(destination, root):
        raise ValueError("Project input STL must remain inside the project directory.")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if not paths_equal(source, destination):
        shutil.copy2(source, destination)
    source_sidecar = source.with_suffix(".geometry.json")
    destination_sidecar = destination.with_suffix(".geometry.json")
    if geometry_manifest is not None:
        shutil.copy2(source_sidecar, destination_sidecar)
    else:
        destination_sidecar.unlink(missing_ok=True)
    region_manifest_path = destination.parent / "region_manifest.json"
    for region_name in ("wall", "inflow", "outflow"):
        (destination.parent / f"{region_name}.stl").unlink(missing_ok=True)
    region_manifest = _copy_region_manifest(source, destination)
    if region_manifest is None:
        region_manifest_path.unlink(missing_ok=True)
    manifest["boundary_regions"] = sorted(region_manifest["regions"]) if region_manifest else None
    manifest["geometry"] = (
        geometry_manifest.get("geometry") if geometry_manifest else None
    )
    manifest["source_stl"] = str(source)
    manifest["template_case_id"] = None
    manifest["mesh_case_id"] = _GENERIC_3D_VALIDATION_CASE
    return _write_manifest(root, manifest)


def save_mesh_project_settings(project: dict, settings: dict) -> dict:
    """Persist editable generation settings without changing the source template."""
    root = Path(project["directory"]).resolve()
    manifest_path = root / "project.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key in ("mesh_size", "mesh_size_min", "mesh_size_max", "boundary_name", "volume_name"):
        if key in settings:
            manifest[key] = settings[key]
    output = normalize_path(settings.get("output_path", project["output_path"]))
    if is_path_within(output, root):
        manifest["output_mesh"] = relative_path_within(output, root).as_posix()
    else:
        manifest["output_mesh"] = str(output)
    return _write_manifest(root, manifest)


__all__ = [
    "create_mesh_project_from_stl",
    "create_mesh_project_from_template",
    "import_mesh_project_directory",
    "load_mesh_project",
    "mesh_project_templates",
    "replace_project_stl",
    "save_mesh_project_settings",
]
