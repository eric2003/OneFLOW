"""Dedicated end-to-end mesh-to-CFD smoke validation."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from grid.stl.stl_to_gmsh import inspect_msh4
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


def mesh_case_choices(repository_root: Path) -> tuple[tuple[str, str], ...]:
    """List native STL-to-Gmsh Cases from the repository Case Catalog."""
    entries = _catalog(Path(repository_root).resolve()).get("cases", [])
    choices = tuple(
        (str(item["id"]), str(item.get("title", item["id"])))
        for item in entries
        if item.get("group") == "STL to Gmsh"
        and (item.get("execution") or {}).get("kind", "oneflow") == "oneflow"
    )
    if not choices:
        raise ValueError("No native STL-to-Gmsh Cases are registered.")
    return choices


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
    if mesh.get("group") != "STL to Gmsh" or (mesh.get("execution") or {}).get("kind", "oneflow") != "oneflow":
        raise ValueError(f"Case {mesh_case_id!r} is not a native STL-to-Gmsh Case.")
    template = next((x for x in _CFD_TEMPLATES if x["id"] == cfd_template_id), None)
    if template is None:
        raise ValueError(f"Unknown CFD template: {cfd_template_id!r}")

    mesh_dir = root / mesh["directory"]
    msh = _setting(mesh_dir / "script/grid.txt", "targetGridFileName")
    stem = Path(msh).stem
    cgns, native = f"results/{stem}.cgns", f"results/{stem}.oneflow.grid"
    mesh_id, cgns_id, native_id = f"{mesh_case_id}_msh", f"{mesh_case_id}_cgns", f"{mesh_case_id}_native"
    cfd_config = root / template["directory"] / "script/cfd.txt"
    return [
        {"id": mesh_case_id, "directory": mesh["directory"],
         "outputs": [{"id": mesh_id, "path": msh, "parameter_file": "script/grid.txt",
                      "setting": "targetGridFileName"}],
         "validation": {"kind": "gmsh_msh4", "path": msh}},
        {"id": f"{mesh_case_id}_gmsh_to_cgns", "directory": "grid/gmsh_to_cgns",
         "depends_on": [mesh_case_id],
         "inputs": [{"artifact": mesh_id, "setting": "sourceGridFileName",
                     "parameter_file": "script/grid.txt"}],
         "outputs": [{"id": cgns_id, "path": cgns}]},
        {"id": f"{mesh_case_id}_cgns_to_oneflow",
         "directory": "grid/stl_to_gmsh/pipeline/cgns_to_oneflow",
         "depends_on": [f"{mesh_case_id}_gmsh_to_cgns"],
         "inputs": [{"artifact": cgns_id, "setting": "sourceGridFileName",
                     "parameter_file": "script/grid.txt"}],
         "outputs": [{"id": native_id, "path": native}]},
        {"id": f"{mesh_case_id}_{cfd_template_id}_solve",
         "directory": template["directory"],
         "depends_on": [f"{mesh_case_id}_cgns_to_oneflow"],
         "inputs": [{"artifact": native_id, "setting": "gridFileName",
                     "parameter_file": "script/cfd.txt"}],
         "outputs": [
             {"id": f"{mesh_case_id}_residual", "path": _setting(cfd_config, "resFile")},
             {"id": f"{mesh_case_id}_restart", "path": _setting(cfd_config, "flowFile")},
             {"id": f"{mesh_case_id}_visual", "path": _setting(cfd_config, "visualFile")},
         ]},
    ]


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
    return UserDataPaths.discover().grid_layout_test / "grid" / "mesh-validation" / relative_case


def reset_mesh_validation_workspace(path: Path) -> Path:
    """Clear and recreate one Case workspace without touching sibling Cases."""
    workspace = _external(path)
    if workspace.is_symlink():
        raise ValueError(f"Mesh validation workspace must not be a symbolic link: {workspace}")
    if workspace.exists():
        shutil.rmtree(workspace)
    workspace.mkdir(parents=True, exist_ok=False)
    return workspace


def run_mesh_validation(
    oneflow: Path, gmsh: str, mesh_case_id: str,
    cfd_template_id: str = "laminar_flat_plate", max_steps: int = 2,
    stage_root: Path | None = None,
) -> int:
    """Stage and execute the selected mesh producer, converters, and CFD Case."""
    executable = oneflow.expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable}")
    gmsh_path = Path(gmsh).expanduser()
    if not gmsh_path.is_file():
        found = shutil.which(gmsh)
        if not found:
            raise FileNotFoundError(f"Gmsh executable not found: {gmsh}")
        gmsh_path = Path(found).resolve()
    if max_steps < 1:
        raise ValueError("CFD smoke-test steps must be at least one.")

    cases = build_mesh_validation_cases(REPOSITORY_ROOT, mesh_case_id, cfd_template_id)
    plan = build_case_graph_plan(cases, REPOSITORY_ROOT)
    run_root = _external(
        stage_root if stage_root is not None
        else mesh_validation_workspace(REPOSITORY_ROOT, mesh_case_id)
    )
    run_root = reset_mesh_validation_workspace(run_root)
    run_id = run_root.name
    os.environ["PATH"] = str(gmsh_path.parent) + os.pathsep + os.environ.get("PATH", "")
    records = []

    def execute(step):
        if step.case_id.endswith("_solve"):
            apply_runtime_overrides(step.directory / "script/cfd.txt", {
                "maxSteps": max_steps, "nFieldSave": 1, "nVisualSave": 1,
                "nResSave": 1, "idump": 0,
            })
        command = [str(executable), "r", str(step.directory)]
        output: list[str] = []
        def emit(line: str) -> None:
            output.append(line)
            print(f"[{step.case_id}] {line}", flush=True)
        try:
            result = run_process(command, cwd=REPOSITORY_ROOT, on_output=emit,
                                 is_cancelled=lambda: False)
            status, message = result.status, result.message
        except OSError as error:
            status, message = "failed", str(error)

        validation = None
        if status == "succeeded" and step.case_id == mesh_case_id:
            try:
                rel = step.outputs[0].relative_path.as_posix()
                mesh_info = inspect_msh4(step.directory / rel)
                validation = {
                    "kind": "gmsh_msh4", "path": rel,
                    "nodes": mesh_info["nodes"], "elements": mesh_info["elements"],
                    "physical_groups": mesh_info["physical_groups"],
                }
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
        records.append({"case_id": step.case_id, "status": status, "message": message})
        print(f"[{step.case_id}] {status}: {message}", flush=True)
        return status, message

    result = run_case_graph(plan, execute, run_root=run_root)
    report = {
        "schema": "oneflow.mesh-to-cfd-validation", "schema_version": 1,
        "run_id": run_id, "status": result.status, "mesh_case_id": mesh_case_id,
        "cfd_template_id": cfd_template_id, "max_steps": max_steps,
        "run_directory": str(run_root), "steps": records,
        "note": "A square mesh with the flat-plate template is a solver smoke test, not physical flat-plate validation.",
        "finished_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    report_path = run_root / "mesh-validation-report.json"
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                           encoding="utf-8", newline="\n")
    print(f"Mesh-to-CFD report: {report_path}", flush=True)
    print(f"Mesh-to-CFD status: {result.status}", flush=True)
    return 0 if result.status == "succeeded" else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oneflow-executable", required=True, type=Path)
    parser.add_argument("--gmsh", default="gmsh")
    parser.add_argument("--mesh-case-id", required=True)
    parser.add_argument("--cfd-template-id", default="laminar_flat_plate")
    parser.add_argument("--max-steps", type=int, default=2)
    parser.add_argument("--stage-root", type=Path)
    args = parser.parse_args()
    try:
        return run_mesh_validation(args.oneflow_executable, args.gmsh,
                                   args.mesh_case_id, args.cfd_template_id,
                                   args.max_steps, args.stage_root)
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, json.JSONDecodeError) as error:
        parser.exit(2, f"Mesh-to-CFD Validation failed: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
