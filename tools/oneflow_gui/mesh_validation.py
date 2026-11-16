"""Dedicated end-to-end mesh-to-CFD smoke validation."""

from __future__ import annotations

import argparse
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


def mesh_case_groups(
    repository_root: Path,
) -> tuple[tuple[str, tuple[tuple[str, str], ...]], ...]:
    """Group supported mesh producers by spatial dimension for GUI navigation."""
    entries = _catalog(Path(repository_root).resolve()).get("cases", [])
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
        raise ValueError("No supported STL-to-Gmsh Cases are registered.")
    return result


def mesh_case_choices(repository_root: Path) -> tuple[tuple[str, str], ...]:
    """List native and Python STL-to-Gmsh producers in catalog order."""
    entries = _catalog(Path(repository_root).resolve()).get("cases", [])
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


def existing_mesh_case_paths(repository_root: Path, mesh_case_id: str) -> tuple[Path, Path] | None:
    """Resolve a catalog-declared generated mesh below the user's mesh workspace."""
    catalog = _catalog(Path(repository_root).resolve())
    mesh = next((item for item in catalog["cases"] if item.get("id") == mesh_case_id), None)
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
    mesh = next((x for x in catalog["cases"] if x.get("id") == mesh_case_id), None)
    if mesh is None:
        raise ValueError(f"Unknown mesh-generation Case: {mesh_case_id!r}")
    execution = mesh.get("execution") or {}
    implementation = execution.get("kind", "oneflow")
    if mesh.get("group") != "STL to Gmsh" or implementation not in {"oneflow", "python", "existing_mesh"}:
        raise ValueError(f"Case {mesh_case_id!r} is not a supported Mesh-to-CFD validation entry.")
    validation_scope = mesh.get("mesh_validation", True)
    if validation_scope is False:
        raise ValueError(
            f"Case {mesh_case_id!r} is not enabled for Mesh-to-CFD validation."
        )
    if validation_scope not in (True, "mesh_handoff", "cfd_smoke"):
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
    if validation_scope == "mesh_handoff":
        # Handoff-only entries do not have a compatible, fully specified CFD
        # consumer contract. Validate MSH, CGNS conversion and OneFLOW ingestion.
        return steps[:3]
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
        (item for item in _catalog(root)["cases"] if item.get("id") == mesh_case_id),
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


def _mesh_validation_stage_results(result) -> dict[str, dict]:
    """Map workflow graph outcomes onto shared run-report states."""
    stage_names = {
        "stl_to_gmsh": "mesh_generation",
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
    return [
        sys.executable, str(script_path), "--case-directory",
        str(case_directory), "--gmsh", str(gmsh_executable),
    ]


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
        (item for item in catalog["cases"] if item.get("id") == mesh_case_id),
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
        if step.case_id == "oneflow_cfd_solve":
            layout_path = REPOSITORY_ROOT / mesh_entry["directory"] / "script/layout.json"
            if layout_path.is_file():
                layout_data = json.loads(layout_path.read_text(encoding="utf-8"))
                if layout_data.get("mesh_dimension") == 3:
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

        validation = None
        if status == "succeeded" and step.case_id == "stl_to_gmsh":
            try:
                rel = step.outputs[0].relative_path.as_posix()
                mesh_info = inspect_msh4(step.directory / rel)
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
    stages = _mesh_validation_stage_results(result)
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
    report.update({
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
