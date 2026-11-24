"""Dedicated end-to-end mesh-to-CFD smoke validation."""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
import json
import os
import re
import shutil
import sys
import uuid
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from grid.stl_to_gmsh.converter import inspect_msh4
from grid.stl_to_gmsh.run_report import build_run_report, format_run_summary, stage_result, write_run_report
from tools.oneflow_gui.user_data_paths import UserDataPaths
from tools.oneflow_gui.workflow_graph_runner import run_case_graph
from tools.oneflow_gui.workflow_plan import build_case_graph_plan
from tools.oneflow_gui.workflow_process import run_process
from tools.oneflow_gui.workflow_runtime_overrides import apply_runtime_overrides

_CFD_TEMPLATES = (
    {"id": "laminar_flat_plate", "title": "Laminar Flat Plate",
     "directory": "test/cfd/plateuns2dslau2"},
)


def _catalog(root: Path) -> dict:
    value = json.loads((root / "grid/case_catalog.json").read_text(encoding="utf-8"))
    if value.get("schema") != "oneflow.case-catalog":
        raise ValueError("Unsupported OneFLOW Case catalog.")
    return value


def _mesh_validation_entries(catalog: dict) -> list[dict]:
    """Return executable mesh Cases plus explicit validation-only mesh inputs."""
    cases = catalog.get("cases", [])
    inputs = catalog.get("mesh_validation_inputs", [])
    if not isinstance(cases, list) or not isinstance(inputs, list):
        raise ValueError("Case catalog mesh validation entries must be lists.")
    return [*cases, *inputs]


def mesh_case_groups(
    repository_root: Path,
) -> tuple[tuple[str, tuple[tuple[str, str], ...]], ...]:
    """Group mesh producers and existing-mesh profiles by spatial dimension."""
    entries = _mesh_validation_entries(_catalog(Path(repository_root).resolve()))
    grouped: dict[int | None, list[tuple[str, str]]] = {2: [], 3: [], None: []}
    for item in entries:
        if item.get("group") != "STL to Gmsh":
            continue
        if item.get("mesh_validation", True) is False:
            continue
        if (item.get("execution") or {}).get("kind", "oneflow") not in {"oneflow", "python", "existing_mesh"}:
            continue
        dimension = item.get("dimension")
        if dimension not in (2, 3):
            dimension = None
        grouped[dimension].append(
            (str(item["id"]), str(item.get("title", item["id"])))
        )
    labels = ((2, "2D"), (3, "3D"), (None, "Unclassified"))
    result = tuple(
        (label, tuple(grouped[dimension]))
        for dimension, label in labels
        if grouped[dimension]
    )
    if not result:
        raise ValueError("No supported STL-to-Gmsh or existing-mesh entries are registered.")
    return result



def mesh_case_group_categories(
    repository_root: Path,
) -> tuple[tuple[str, tuple[tuple[str, tuple[tuple[str, str], ...]], ...]], ...]:
    """Separate generated STL-to-Gmsh Cases from pre-existing MSH inputs."""
    catalog = _catalog(Path(repository_root).resolve())
    entries = {
        str(item["id"]): item
        for item in _mesh_validation_entries(catalog)
        if isinstance(item, dict) and "id" in item
    }
    result = []
    for dimension, choices in mesh_case_groups(repository_root):
        generated = []
        existing = []
        for choice in choices:
            entry = entries.get(choice[0], {})
            kind = (entry.get("execution") or {}).get("kind", "oneflow")
            (existing if kind == "existing_mesh" else generated).append(choice)
        categories = []
        if generated:
            categories.append(("STL to Gmsh", tuple(generated)))
        if existing:
            categories.append(("Existing Gmsh Mesh", tuple(existing)))
        result.append((dimension, tuple(categories)))
    return tuple(result)

def mesh_case_choices(repository_root: Path) -> tuple[tuple[str, str], ...]:
    """List supported STL-to-Gmsh producers and existing-mesh inputs."""
    entries = _mesh_validation_entries(_catalog(Path(repository_root).resolve()))
    choices = tuple(
        (str(item["id"]), str(item.get("title", item["id"])))
        for item in entries
        if item.get("group") == "STL to Gmsh"
        and item.get("mesh_validation", True) is not False
        and (item.get("execution") or {}).get("kind", "oneflow") in {"oneflow", "python", "existing_mesh"}
    )
    if not choices:
        raise ValueError("No supported STL-to-Gmsh Cases are registered.")
    return choices


def existing_mesh_case_paths(
    repository_root: Path, mesh_case_id: str,
) -> tuple[Path, Path] | None:
    """Resolve a catalog-declared generated mesh below the user's mesh workspace."""
    catalog = _catalog(Path(repository_root).resolve())
    mesh = next((item for item in _mesh_validation_entries(catalog) if item.get("id") == mesh_case_id), None)
    if mesh is None:
        raise ValueError(f"Unknown mesh-generation Case: {mesh_case_id!r}")
    execution = mesh.get("execution") or {}
    if execution.get("kind") != "existing_mesh":
        return None

    project_name = execution.get("project_directory")
    mesh_value = execution.get("mesh_file")
    if not isinstance(project_name, str) or not isinstance(mesh_value, str):
        raise ValueError(f"Existing-mesh Case {mesh_case_id!r} must declare its project and mesh path.")
    project_relative = Path(project_name)
    mesh_relative = Path(mesh_value)
    if (
        project_relative.is_absolute() or not project_relative.parts
        or ".." in project_relative.parts
        or mesh_relative.is_absolute() or not mesh_relative.parts
        or ".." in mesh_relative.parts
        or mesh_relative.suffix.lower() != ".msh"
    ):
        raise ValueError(f"Existing-mesh paths for {mesh_case_id!r} must be safe relative .msh paths.")

    workspace = UserDataPaths.discover().mesh_generation.resolve()
    project_directory = (workspace / project_relative).resolve()
    mesh_path = (project_directory / mesh_relative).resolve()
    try:
        project_directory.relative_to(workspace)
        mesh_path.relative_to(project_directory)
    except ValueError as error:
        raise ValueError(f"Existing-mesh paths for {mesh_case_id!r} escape the mesh workspace.") from error
    return mesh_path, project_directory


def cfd_template_choices() -> tuple[tuple[str, str], ...]:
    """List available CFD consumer templates."""
    return tuple((item["id"], item["title"]) for item in _CFD_TEMPLATES)


def _setting(path: Path, name: str) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf'^\s*string\s+{re.escape(name)}\s*=\s*"([^"]+)"\s*;',
        text, re.MULTILINE,
    )
    if match is None:
        raise ValueError(f"Missing string setting {name!r} in {path}.")
    value = Path(match.group(1))
    if value.is_absolute() or ".." in value.parts:
        raise ValueError(f"Setting {name!r} must be a Case-relative path.")
    return value.as_posix()


def build_mesh_validation_cases(
    repository_root: Path, mesh_case_id: str, cfd_template_id: str,
) -> list[dict]:
    """Build the current STL -> Gmsh -> CGNS -> OneFLOW -> CFD route."""
    root = Path(repository_root).resolve()
    catalog = _catalog(root)
    mesh = next((x for x in _mesh_validation_entries(catalog) if x.get("id") == mesh_case_id), None)
    if mesh is None:
        raise ValueError(f"Unknown mesh-generation Case: {mesh_case_id!r}")
    execution = mesh.get("execution") or {}
    implementation = execution.get("kind", "oneflow")
    if (
        mesh.get("group") != "STL to Gmsh"
        or implementation not in {"oneflow", "python", "existing_mesh"}
    ):
        raise ValueError(
            f"Case {mesh_case_id!r} is not a supported Mesh-to-CFD validation entry."
        )
    validation_scope = mesh.get("mesh_validation", True)
    if validation_scope is False:
        raise ValueError(
            f"Case {mesh_case_id!r} is not enabled for Mesh-to-CFD validation."
        )
    if validation_scope not in (True, "cfd_smoke"):
        raise ValueError(
            f"Unsupported Mesh-to-CFD validation scope for {mesh_case_id!r}: "
            f"{validation_scope!r}"
        )
    template = next((x for x in _CFD_TEMPLATES if x["id"] == cfd_template_id), None)
    if template is None:
        raise ValueError(f"Unknown CFD template: {cfd_template_id!r}")

    mesh_dir = root / mesh["directory"]
    if implementation == "oneflow":
        msh = _setting(mesh_dir / "script/grid.txt", "targetGridFileName")
    else:
        expected_path = mesh_dir / "expected.json"
        expected = json.loads(expected_path.read_text(encoding="utf-8"))
        mesh_name = expected.get("artifact", {}).get("mesh_file")
        if (
            not isinstance(mesh_name, str)
            or not mesh_name.endswith(".msh")
            or Path(mesh_name).name != mesh_name
            or "\\" in mesh_name
        ):
            raise ValueError(
                f"Mesh Case {mesh_case_id!r} must declare artifact.mesh_file "
                "as a simple .msh filename."
            )
        msh = f"results/{mesh_name}"
    stem = Path(msh).stem
    cgns, native = f"results/{stem}.cgns", f"results/{stem}.oneflow.grid"
    bc_region_map = f"results/{stem}_bcRegionMap.txt"
    grid_info = f"results/{stem}.info.json"
    mesh_id, cgns_id, native_id = "gmsh_mesh", "cgns_mesh", "oneflow_grid"
    cfd_config = root / template["directory"] / "script/cfd.txt"
    steps = [
        {"id": "stl_to_gmsh", "directory": mesh["directory"],
         "runtime_directory": "stl_to_gmsh",
         "outputs": (
             [{"id": mesh_id, "path": msh, "parameter_file": "script/grid.txt",
               "setting": "targetGridFileName"}]
             if implementation == "oneflow"
             else [{"id": mesh_id, "path": msh}]
         ),
         "validation": {"kind": "gmsh_msh4", "path": msh}},
        {"id": "gmsh_to_cgns", "directory": "grid/conversion/gmsh_to_cgns",
         "runtime_directory": "gmsh_to_cgns",
         "depends_on": ["stl_to_gmsh"],
         "inputs": [{"artifact": mesh_id, "setting": "sourceGridFileName",
                     "parameter_file": "script/grid.txt"}],
         "outputs": [{"id": cgns_id, "path": cgns,
                      "parameter_file": "script/grid.txt",
                      "setting": "targetGridFileName"}]},
        {"id": "cgns_to_oneflow",
         "directory": "grid/conversion/cgns_to_oneflow",
         "runtime_directory": "cgns_to_oneflow",
         "depends_on": ["gmsh_to_cgns"],
         "inputs": [{"artifact": cgns_id, "setting": "sourceGridFileName",
                     "parameter_file": "script/grid.txt"}],
         "outputs": [
             {"id": native_id, "path": native,
              "parameter_file": "script/grid.txt",
              "setting": "targetGridFileName"},
             {"id": "bc_region_map", "path": bc_region_map,
              "parameter_file": "script/grid.txt",
              "setting": "bcRegionMapFile"},
             {"id": "grid_info", "path": grid_info,
              "parameter_file": "script/grid.txt",
              "setting": "gridInfoFileName"},
         ]},
        {"id": "oneflow_cfd_solve",
         "directory": template["directory"],
         "runtime_directory": "oneflow_cfd_solve",
         "depends_on": ["cgns_to_oneflow"],
         "inputs": [
             {"artifact": native_id, "setting": "gridFileName",
              "parameter_file": "script/cfd.txt"},
             {"artifact": "bc_region_map", "setting": "bcRegionMapFile",
              "parameter_file": "script/cfd.txt"},
         ],
         "outputs": [
             {"id": "residual", "path": _setting(cfd_config, "resFile")},
             {"id": "restart", "path": _setting(cfd_config, "flowFile")},
             {"id": "visual", "path": _setting(cfd_config, "visualFile")},
         ]},
    ]
    return steps


def _external(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(REPOSITORY_ROOT)
    except ValueError:
        return resolved
    raise ValueError(f"Run data must stay outside the source repository: {resolved}")


def mesh_validation_workspace(repository_root: Path, mesh_case_id: str) -> Path:
    """Return the stable user-data directory for one mesh validation Case."""
    root = Path(repository_root).resolve()
    entry = next(
        (item for item in _mesh_validation_entries(_catalog(root)) if item.get("id") == mesh_case_id),
        None,
    )
    if entry is None:
        raise ValueError(f"Unknown mesh-generation Case: {mesh_case_id!r}")
    case_directory = Path(entry["directory"])
    if (
        case_directory.is_absolute()
        or ".." in case_directory.parts
        or len(case_directory.parts) < 2
        or case_directory.parts[0] != "grid"
    ):
        raise ValueError(
            f"Mesh validation Case must be located below grid/: {case_directory}"
        )
    relative_case = Path(*case_directory.parts[1:])
    return UserDataPaths.discover().mesh_validation / relative_case


def project_mesh_validation_workspace(project_directory: Path) -> Path:
    """Return the stable validation directory belonging to one generated mesh project."""
    project_root = _external(project_directory)
    return project_root / "validation"


def reset_mesh_validation_workspace(path: Path) -> Path:
    """Clear and recreate one validation workspace without touching its mesh project."""
    workspace = _external(path)
    if workspace.is_symlink():
        raise ValueError(f"Mesh validation workspace must not be a symbolic link: {workspace}")
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True, exist_ok=False)
    return workspace


def create_mesh_validation_run_workspace(path: Path) -> Path:
    """Clear and recreate the validation directory before each new run."""
    return reset_mesh_validation_workspace(path)


def _graph_step_records(result) -> list[dict[str, str]]:
    """Serialize final graph outcomes, including artifact checks and blocked steps."""
    return [
        {"case_id": item.case_id, "status": item.status, "message": item.message}
        for item in result.cases
    ]


def run_report_stage_rows(report: dict) -> tuple[tuple[str, str, str], ...]:
    """Return display-ready stage rows from a shared run report."""
    stages = report.get("stages", {})
    if not isinstance(stages, dict):
        raise ValueError("Run report stages must be an object.")
    return tuple(
        (
            str(name),
            str(result.get("status", "UNKNOWN")),
            str(result.get("message", "")),
        )
        for name, result in stages.items()
        if isinstance(result, dict)
    )


def geometry_lineage_rows(report: dict) -> tuple[tuple[str, str, str], ...]:
    """Return readable geometry-source and mesh-derivation rows from a run report."""
    lineage = report.get("geometry_lineage") or {}
    if not isinstance(lineage, dict):
        return ()
    provenance = lineage.get("provenance") or {}
    kind = str(provenance.get("kind", "unknown")).upper()
    notes = str(provenance.get("notes", ""))
    rows = [("Geometry provenance", kind, notes or "No verified source record is available.")]
    asset_id = lineage.get("asset_id")
    if asset_id:
        rows.append(("Geometry asset", "IDENTIFIED", str(asset_id)))
    if lineage.get("manifest_status"):
        rows.append(("Geometry manifest", "UNTRUSTED", str(lineage["manifest_status"])))
    source = lineage.get("source_representation") or {}
    if source:
        rows.append(("Source surface", "RECORDED",
                     f"{source.get('format', 'unknown')} | {source.get('path', 'unknown')} | SHA-256: {source.get('sha256', 'unknown')}"))
    generator = lineage.get("generator")
    if isinstance(generator, dict):
        if generator:
            name = str(generator.get("name", "unknown"))
            version = str(generator.get("version", "")).strip()
            source_sha256 = generator.get("source_sha256")
            details = " ".join(part for part in (name, version) if part)
            if source_sha256:
                details += f" | source SHA-256: {source_sha256}"
            rows.append(("Geometry generator", "RECORDED", details))
    elif isinstance(generator, str) and generator.strip():
        # Older and lightweight manifests may store the generator as a name.
        rows.append(("Geometry generator", "RECORDED", generator.strip()))
    elif generator is not None:
        rows.append(("Geometry generator", "UNRECOGNIZED", str(generator)))
    mesh = lineage.get("derived_mesh") or {}
    if mesh:
        rows.append(("Derived volume mesh", "RECORDED",
                     f"{mesh.get('format', 'unknown')} | {mesh.get('path', 'unknown')} | SHA-256: {mesh.get('sha256', 'unknown')}"))
    for artifact in lineage.get("downstream_artifacts", []):
        rows.append((
            str(artifact.get("role", "Derived artifact")),
            "RECORDED",
            f"{artifact.get('format', 'unknown')} | {artifact.get('path', 'unknown')} | SHA-256: {artifact.get('sha256', 'unknown')}",
        ))
    return tuple(rows)


def _geometry_lineage(
    mesh_entry: dict, source_stl: Path | None, mesh_file: Path | None,
    repository_root: Path, downstream_artifacts: list[dict] | None = None,
) -> dict:
    """Capture only verifiable source metadata and hashes for this validation run."""
    provenance = mesh_entry.get("geometry_provenance") or {
        "kind": "unknown",
        "notes": "No verified geometry provenance is registered for this source Case.",
    }
    lineage = {"provenance": provenance}
    manifest = None
    if source_stl is not None and source_stl.is_file():
        source_stl = source_stl.resolve()
        digest = hashlib.sha256(source_stl.read_bytes()).hexdigest()
        try:
            source_label = source_stl.relative_to(repository_root.resolve()).as_posix()
        except ValueError:
            source_label = str(source_stl)
        lineage["source_representation"] = {
            "format": "STL", "path": source_label, "sha256": digest,
        }
        manifest_path = source_stl.with_suffix(".geometry.json")
        if manifest_path.is_file():
            try:
                candidate = json.loads(manifest_path.read_text(encoding="utf-8"))
                if (isinstance(candidate, dict)
                        and candidate.get("schema") == "oneflow.stl-geometry-manifest"
                        and candidate.get("stl_sha256") == digest):
                    manifest = candidate
                else:
                    lineage["manifest_status"] = "invalid_or_checksum_mismatch"
            except (OSError, json.JSONDecodeError):
                lineage["manifest_status"] = "unreadable"
    if manifest is not None:
        lineage["asset_id"] = manifest.get("asset_id")
        lineage["definition_sha256"] = manifest.get("definition_sha256")
        lineage["generator"] = manifest.get("generator", {})
        lineage["provenance"] = {
            "kind": "generated",
            "notes": "Verified companion geometry manifest matches the source STL checksum.",
        }
    if mesh_file is not None and mesh_file.is_file():
        mesh_file = mesh_file.resolve()
        try:
            mesh_label = mesh_file.relative_to(repository_root.resolve()).as_posix()
        except ValueError:
            mesh_label = str(mesh_file)
        lineage["derived_mesh"] = {
            "format": "Gmsh MSH", "path": mesh_label,
            "sha256": hashlib.sha256(mesh_file.read_bytes()).hexdigest(),
            "parent": "source_representation",
        }
    downstream = []
    for item in downstream_artifacts or []:
        path = Path(item["path"]).expanduser().resolve()
        if not path.is_file():
            continue
        try:
            label = path.relative_to(repository_root.resolve()).as_posix()
        except ValueError:
            label = str(path)
        downstream.append({
            "role": str(item.get("role", "derived artifact")),
            "format": str(item.get("format", "unknown")),
            "path": label,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "parent": item.get("parent", "derived_mesh"),
        })
    if downstream:
        lineage["downstream_artifacts"] = downstream
    return lineage




def _geometry_artifact_lineage_path(artifact_path: Path) -> Path:
    """Return the companion lineage path without changing the artifact itself."""
    return artifact_path.with_name(artifact_path.name + ".geometry-lineage.json")


def _write_geometry_artifact_lineage(
    artifact_path: Path, role: str, geometry_lineage: dict,
    parent_artifact: dict | None = None,
) -> Path:
    """Persist a checksum-bound lineage record beside a generated artifact."""
    artifact = Path(artifact_path).expanduser().resolve()
    if not artifact.is_file():
        raise FileNotFoundError(f"Cannot attach geometry lineage to missing artifact: {artifact}")
    record = {
        "schema": "oneflow.geometry-artifact-lineage",
        "version": 1,
        "artifact": {
            "role": str(role),
            "file_name": artifact.name,
            "format": artifact.suffix.lstrip(".").upper() or "unknown",
            "sha256": hashlib.sha256(artifact.read_bytes()).hexdigest(),
        },
        "geometry_lineage": geometry_lineage,
        "parent_artifact": parent_artifact,
    }
    sidecar = _geometry_artifact_lineage_path(artifact)
    sidecar.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n",
    )
    return sidecar


def _read_geometry_artifact_lineage(artifact_path: Path) -> dict:
    """Read lineage only when the companion record matches the artifact bytes."""
    artifact = Path(artifact_path).expanduser().resolve()
    sidecar = _geometry_artifact_lineage_path(artifact)
    record = json.loads(sidecar.read_text(encoding="utf-8"))
    if not isinstance(record, dict) or record.get("schema") != "oneflow.geometry-artifact-lineage":
        raise ValueError(f"Unsupported geometry lineage sidecar: {sidecar}")
    metadata = record.get("artifact")
    if not isinstance(metadata, dict):
        raise ValueError(f"Geometry lineage sidecar has no artifact record: {sidecar}")
    actual = hashlib.sha256(artifact.read_bytes()).hexdigest()
    if metadata.get("file_name") != artifact.name:
        raise ValueError(f"Geometry lineage filename mismatch for artifact: {artifact}")
    if metadata.get("sha256") != actual:
        raise ValueError(f"Geometry lineage checksum mismatch for artifact: {artifact}")
    lineage = record.get("geometry_lineage")
    if not isinstance(lineage, dict):
        raise ValueError(f"Geometry lineage sidecar has no geometry lineage: {sidecar}")
    return record


def _mesh_validation_stage_results(result, *, existing_mesh: bool = False) -> dict[str, dict]:
    """Map workflow graph outcomes onto shared run-report states."""
    stage_names = {
        "stl_to_gmsh": "existing_msh_input" if existing_mesh else "mesh_generation",
        "gmsh_to_cgns": "gmsh_to_cgns",
        "cgns_to_oneflow": "cgns_to_oneflow",
        "oneflow_cfd_solve": "cfd_smoke",
    }
    stages: dict[str, dict] = {}
    for case in result.cases:
        name = stage_names.get(case.case_id)
        if name is None:
            continue
        status = {
            "succeeded": "PASS", "failed": "FAIL",
            "skipped": "NOT_RUN", "cancelled": "NOT_RUN",
        }.get(case.status)
        if status is None:
            raise ValueError(f"Unsupported workflow result for {case.case_id}: {case.status!r}")
        stages[name] = stage_result(status, case.message)
    for name in stage_names.values():
        stages.setdefault(name, stage_result("NOT_RUN", "Stage was not reached."))
    stages["mesh_quality"] = stage_result(
        "NOT_CHECKED", "No mesh-quality acceptance thresholds are applied yet.")
    stages["cfd_accuracy"] = stage_result(
        "NOT_CHECKED", "A runtime smoke test does not establish physical accuracy.")
    return stages


def _mesh_case_command(
    repository_root: Path,
    mesh_entry: dict,
    oneflow_executable: Path,
    gmsh_executable: Path,
    case_directory: Path,
) -> list[str]:
    """Build the selected mesh producer command without coupling to its language."""
    execution = mesh_entry.get("execution") or {}
    kind = execution.get("kind", "oneflow")
    if kind == "oneflow":
        return [str(oneflow_executable), "r", str(case_directory)]
    if kind != "python":
        raise ValueError(f"Unsupported mesh producer implementation: {kind!r}")
    script_value = execution.get("script")
    if not isinstance(script_value, str) or not script_value.strip():
        raise ValueError("Python mesh producer must declare execution.script.")
    root = repository_root.resolve()
    script_path = (root / script_value).resolve()
    try:
        script_path.relative_to(root)
    except ValueError as error:
        raise ValueError("Python mesh producer script must stay inside the repository.") from error
    command = [
        sys.executable, str(script_path), "--case-directory",
        str(case_directory), "--gmsh", str(gmsh_executable),
    ]
    output_argument = execution.get("output_directory_argument")
    if output_argument:
        if not isinstance(output_argument, str) or not output_argument.startswith("--"):
            raise ValueError("Python mesh producer output_directory_argument must be a CLI option.")
        command.extend([output_argument, str(case_directory)])
    return command


def _project_required_physical_groups(project_data: dict) -> list[dict]:
    """Translate checksum-verified project boundary-region names into MSH requirements."""
    regions = project_data.get("boundary_regions")
    if regions is None:
        return []
    if (
        not isinstance(regions, list)
        or not regions
        or any(not isinstance(name, str) or not name.strip() for name in regions)
        or len(set(regions)) != len(regions)
    ):
        raise ValueError("Project boundary_regions must be a non-empty list of unique names.")
    if any(name not in {"wall", "inflow", "outflow"} for name in regions):
        raise ValueError("Project boundary_regions contains an unsupported STL region name.")
    return [
        *({"name": name, "dimension": 2} for name in regions),
        {"name": "fluid_domain", "dimension": 3},
    ]


def _validate_required_physical_groups(mesh_info: dict, mesh_entry: dict) -> None:
    """Reject an MSH that omits declared boundary or volume Physical Groups."""
    required = mesh_entry.get("required_physical_groups", [])
    if not required:
        return
    available = {
        (int(group.get("dimension", -1)), str(group.get("name", "")))
        for group in mesh_info.get("physical_groups", [])
        if isinstance(group, dict)
    }
    missing = []
    for group in required:
        if (
            not isinstance(group, dict)
            or not isinstance(group.get("name"), str)
            or group.get("dimension") not in (1, 2, 3)
        ):
            raise ValueError("Case required_physical_groups contains an invalid declaration.")
        key = (int(group["dimension"]), group["name"])
        if key not in available:
            missing.append(f"{group['name']} (dimension {group['dimension']})")
    if missing:
        raise ValueError(
            "MSH is missing required Physical Groups: " + ", ".join(missing)
        )


def _copy_existing_mesh(source: Path, destination: Path) -> None:
    """Copy a user-selected MSH file into an isolated workflow stage."""
    source = source.expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Existing Gmsh mesh not found: {source}")
    if source.suffix.lower() != ".msh":
        raise ValueError(f"Existing mesh must use the .msh extension: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source == destination.resolve():
        raise ValueError("Existing mesh source and staged destination must differ.")
    shutil.copy2(source, destination)


def _downstream_lineage_artifacts(steps: list[dict], run_root: Path) -> list[dict]:
    """Resolve CGNS and native-grid outputs for the final geometry lineage report."""
    artifacts = []
    for step in steps:
        if step.get("id") not in {"gmsh_to_cgns", "cgns_to_oneflow"}:
            continue
        for output in step.get("outputs", []):
            artifact_id = output.get("id")
            if artifact_id not in {"cgns_mesh", "oneflow_grid"}:
                continue
            artifacts.append({
                "role": "CGNS volume mesh" if artifact_id == "cgns_mesh" else "OneFLOW native grid",
                "format": "CGNS" if artifact_id == "cgns_mesh" else "OneFLOW native grid",
                "path": run_root / step["runtime_directory"] / output["path"],
                "parent": "derived_mesh" if artifact_id == "cgns_mesh" else "cgns_mesh",
            })
    return artifacts


def run_mesh_validation(
    oneflow: Path, gmsh: str, mesh_case_id: str,
    cfd_template_id: str = "laminar_flat_plate", max_steps: int = 2,
    stage_root: Path | None = None, mesh_file: Path | None = None,
    project_directory: Path | None = None,
) -> int:
    """Run a mesh producer or validate an existing MSH, then convert and consume it."""
    executable = oneflow.expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable}")
    source_mesh = mesh_file.expanduser().resolve() if mesh_file is not None else None
    if max_steps < 1:
        raise ValueError("CFD smoke-test steps must be at least one.")

    catalog = _catalog(REPOSITORY_ROOT)
    mesh_entry = next(
        (item for item in _mesh_validation_entries(catalog) if item.get("id") == mesh_case_id),
        None,
    )
    if mesh_entry is None:
        raise ValueError(f"Unknown mesh-generation Case: {mesh_case_id!r}")
    mesh_execution = mesh_entry.get("execution") or {}
    if source_mesh is None and mesh_execution.get("kind") == "existing_mesh":
        existing_mesh_paths = existing_mesh_case_paths(REPOSITORY_ROOT, mesh_case_id)
        if existing_mesh_paths is None:
            raise ValueError(f"Existing-mesh Case {mesh_case_id!r} has no mesh path.")
        source_mesh, default_project_directory = existing_mesh_paths
        if project_directory is None:
            project_directory = default_project_directory

    source_stl = None
    if project_directory is not None:
        project_root = Path(project_directory).expanduser().resolve()
        project_manifest = project_root / "project.json"
        if project_manifest.is_file():
            project_data = json.loads(project_manifest.read_text(encoding="utf-8"))
            project_groups = _project_required_physical_groups(project_data)
            if project_groups:
                mesh_entry = {**mesh_entry, "required_physical_groups": project_groups}
            relative_stl = Path(str(project_data.get("input_stl", "")))
            if not relative_stl.is_absolute() and ".." not in relative_stl.parts:
                candidate = (project_root / relative_stl).resolve()
                try:
                    candidate.relative_to(project_root)
                except ValueError:
                    candidate = None
                if candidate is not None and candidate.is_file():
                    source_stl = candidate
    if source_stl is None and source_mesh is None and project_directory is None:
        candidate = REPOSITORY_ROOT / mesh_entry["directory"] / "input" / "model.stl"
        if candidate.is_file():
            source_stl = candidate

    gmsh_path = Path(gmsh).expanduser()
    if source_mesh is None and not gmsh_path.is_file():
        found = shutil.which(gmsh)
        if not found:
            raise FileNotFoundError(f"Gmsh executable not found: {gmsh}")
        gmsh_path = Path(found).resolve()
    mesh_implementation = (
        "existing_mesh" if source_mesh is not None
        else mesh_execution.get("kind", "oneflow")
    )
    if source_mesh is not None:
        if not source_mesh.is_file():
            raise FileNotFoundError(f"Existing Gmsh mesh not found: {source_mesh}")
        if source_mesh.suffix.lower() != ".msh":
            raise ValueError(f"Existing mesh must use the .msh extension: {source_mesh}")
        inspect_msh4(source_mesh)
    cases = build_mesh_validation_cases(REPOSITORY_ROOT, mesh_case_id, cfd_template_id)
    plan = build_case_graph_plan(cases, REPOSITORY_ROOT)
    if stage_root is not None:
        selected_stage_root = stage_root
    elif source_mesh is not None:
        selected_stage_root = project_mesh_validation_workspace(
            project_directory or source_mesh.parent
        )
    else:
        selected_stage_root = mesh_validation_workspace(REPOSITORY_ROOT, mesh_case_id)
    run_root = _external(selected_stage_root)
    run_root = reset_mesh_validation_workspace(run_root)
    run_id = uuid.uuid4().hex
    mesh_evidence: dict = {}
    if gmsh_path.is_file():
        os.environ["PATH"] = str(gmsh_path.parent) + os.pathsep + os.environ.get("PATH", "")
    def execute(step):
        inherited_cgns_lineage = None
        if step.case_id == "oneflow_cfd_solve":
            # Dimension is part of the mesh Case contract. Do not infer it from
            # an optional layout.json: existing-MSH projects need the same 3D
            # solver setup as generated sphere and spheroid Cases.
            if mesh_entry.get("dimension") == 3:
                apply_runtime_overrides(
                    step.directory / "script/dimension.txt", {"dimension": 3}
                )
            apply_runtime_overrides(step.directory / "script/cfd.txt", {
                "maxSteps": max_steps, "nFieldSave": 1, "nVisualSave": 1,
                "nResSave": 1, "idump": 0,
            })
        if step.case_id == "stl_to_gmsh" and source_mesh is not None:
            command = ["<copy-existing-mesh>", str(source_mesh)]
            try:
                relative_mesh = step.outputs[0].relative_path
                _copy_existing_mesh(source_mesh, step.directory / relative_mesh)
                status, message = "succeeded", f"Staged existing mesh: {source_mesh}"
            except (OSError, ValueError) as error:
                status, message = "failed", str(error)
        elif step.case_id == "stl_to_gmsh":
            command = _mesh_case_command(
                REPOSITORY_ROOT, mesh_entry, executable, gmsh_path, step.directory
            )
            status, message = "", ""
        else:
            command = [str(executable), "r", str(step.directory)]
            status, message = "", ""
        if step.case_id == "cgns_to_oneflow":
            try:
                cgns_input = step.inputs[0].source_path
                inherited_cgns_lineage = _read_geometry_artifact_lineage(cgns_input)
            except (OSError, ValueError, json.JSONDecodeError, IndexError) as error:
                status, message = "failed", f"CGNS geometry lineage verification failed: {error}"

        output: list[str] = []
        def emit(line: str) -> None:
            output.append(line)
            print(f"[{step.case_id}] {line}", flush=True)
        if not status:
            try:
                result = run_process(command, cwd=step.directory, on_output=emit,
                                     is_cancelled=lambda: False)
                status, message = result.status, result.message
            except OSError as error:
                status, message = "failed", str(error)

        if status == "succeeded" and step.case_id == "gmsh_to_cgns":
            try:
                msh_input = step.inputs[0].source_path
                geometry_lineage = _geometry_lineage(
                    mesh_entry, source_stl, msh_input, REPOSITORY_ROOT
                )
                cgns_output = step.directory / next(
                    item.relative_path for item in step.outputs
                    if item.artifact_id == "cgns_mesh"
                )
                sidecar = _write_geometry_artifact_lineage(
                    cgns_output, "CGNS volume mesh", geometry_lineage,
                    {
                        "role": "Gmsh MSH", "file_name": msh_input.name,
                        "path": msh_input.relative_to(run_root).as_posix(),
                        "sha256": hashlib.sha256(msh_input.read_bytes()).hexdigest(),
                    },
                )
                print(
                    f"[{step.case_id}] geometry lineage: {sidecar.name}",
                    flush=True,
                )
            except (OSError, ValueError, StopIteration) as error:
                status, message = "failed", f"CGNS geometry lineage persistence failed: {error}"

        if status == "succeeded" and step.case_id == "cgns_to_oneflow":
            try:
                cgns_input = step.inputs[0].source_path
                native_output = step.directory / next(
                    item.relative_path for item in step.outputs
                    if item.artifact_id == "oneflow_grid"
                )
                cgns_artifact = inherited_cgns_lineage["artifact"]
                sidecar = _write_geometry_artifact_lineage(
                    native_output, "OneFLOW native grid",
                    inherited_cgns_lineage["geometry_lineage"],
                    {
                        "role": "CGNS volume mesh",
                        "file_name": cgns_input.name,
                        "path": cgns_input.relative_to(run_root).as_posix(),
                        "sha256": cgns_artifact["sha256"],
                    },
                )
                print(
                    f"[{step.case_id}] geometry lineage: {sidecar.name}",
                    flush=True,
                )
            except (OSError, ValueError, StopIteration, TypeError, KeyError) as error:
                status, message = "failed", f"OneFLOW grid geometry lineage persistence failed: {error}"

        validation = None
        if status == "succeeded" and step.case_id == "stl_to_gmsh":
            try:
                rel = step.outputs[0].relative_path.as_posix()
                mesh_info = inspect_msh4(step.directory / rel)
                _validate_required_physical_groups(mesh_info, mesh_entry)
                validation = {
                    "kind": "gmsh_msh4", "path": rel,
                    "implementation": mesh_implementation,
                    "nodes": mesh_info["nodes"], "elements": mesh_info["elements"],
                    "physical_groups": mesh_info["physical_groups"],
                }
                mesh_evidence.update(validation)
            except (OSError, ValueError, KeyError, RuntimeError) as error:
                status, message = "failed", f"Gmsh MSH4 validation failed: {error}"

        results = step.directory / "results"
        results.mkdir(parents=True, exist_ok=True)
        (results / "mesh-validation.log").write_text(
            "$ " + json.dumps(command, ensure_ascii=True) + "\n" + "\n".join(output) + "\n",
            encoding="utf-8", newline="\n",
        )
        if validation is not None:
            (results / "mesh-validation.json").write_text(
                json.dumps(validation, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8", newline="\n",
            )
        print(f"[{step.case_id}] {status}: {message}", flush=True)
        return status, message

    result = run_case_graph(plan, execute, run_root=run_root)
    stages = _mesh_validation_stage_results(
        result, existing_mesh=mesh_implementation == "existing_mesh"
    )
    exit_code = 0 if result.status == "succeeded" else 130 if result.status == "cancelled" else 1
    failed_step = next((item for item in result.cases if item.status == "failed"), None)
    error = None
    if failed_step is not None:
        error = {"type": "WorkflowStageError",
                 "message": f"{failed_step.case_id}: {failed_step.message}"}
    elif result.status == "cancelled":
        error = {"type": "WorkflowCancelled", "message": "The validation run was cancelled."}
    metrics = {key: mesh_evidence[key] for key in ("nodes", "elements") if key in mesh_evidence}
    if "physical_groups" in mesh_evidence:
        metrics["physical_group_count"] = len(mesh_evidence["physical_groups"])

    report_path = run_root / "mesh-validation-report.json"
    report = build_run_report(
        "Mesh-to-CFD Validation", stages, exit_code, metrics=metrics,
        artifacts={"run_directory": str(run_root), "report": str(report_path),
                   **({"source_mesh": str(source_mesh)} if source_mesh is not None else {})},
        error=error,
    )
    staged_mesh = None
    if mesh_evidence.get("path"):
        candidate = run_root / "stl_to_gmsh" / str(mesh_evidence["path"])
        if candidate.is_file():
            staged_mesh = candidate
    downstream_artifacts = _downstream_lineage_artifacts(cases, run_root)
    lineage = _geometry_lineage(
        mesh_entry, source_stl, staged_mesh, REPOSITORY_ROOT, downstream_artifacts
    )
    report.update({
        "geometry_lineage": lineage,
        "legacy_schema": "oneflow.mesh-to-cfd-validation",
        "schema_version": 1, "run_id": run_id, "mesh_case_id": mesh_case_id,
        "mesh_implementation": mesh_implementation,
        "mesh_input_mode": "existing_mesh" if source_mesh is not None else "generated_mesh",
        "source_mesh_file": str(source_mesh) if source_mesh is not None else None,
        "cfd_template_id": cfd_template_id, "max_steps": max_steps,
        "run_directory": str(run_root), "steps": _graph_step_records(result),
        "mesh_summary": mesh_evidence or None,
        "note": "This report records execution and ingestion smoke checks only; it does not establish mesh-quality acceptance or CFD physical accuracy.",
        "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    })
    write_run_report(report_path, report)
    print(format_run_summary(report), flush=True)
    return exit_code


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oneflow-executable", required=True, type=Path)
    parser.add_argument("--gmsh", default="gmsh")
    parser.add_argument("--mesh-case-id", required=True)
    parser.add_argument("--cfd-template-id", default="laminar_flat_plate")
    parser.add_argument("--max-steps", type=int, default=2)
    parser.add_argument("--stage-root", type=Path)
    parser.add_argument("--mesh-file", type=Path,
                        help="Validate an existing MSH4 file without rerunning its producer.")
    parser.add_argument("--project-directory", type=Path,
                        help="Mesh project root; validation stages are written directly below its validation/ directory.")
    args = parser.parse_args()
    try:
        return run_mesh_validation(args.oneflow_executable, args.gmsh,
                                   args.mesh_case_id, args.cfd_template_id,
                                   args.max_steps, args.stage_root,
                                   args.mesh_file, args.project_directory)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, json.JSONDecodeError) as error:
        parser.exit(2, f"Mesh-to-CFD Validation failed: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
