#!/usr/bin/env python3
"""Run a declarative STL -> Gmsh -> CGNS -> OneFLOW -> CFD case graph."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from uuid import uuid4

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from grid.stl.stl_to_gmsh import inspect_msh4
from tools.oneflow_gui.run_environment_evidence import RunEnvironmentEvidence
from tools.oneflow_gui.workflow_graph_runner import run_case_graph
from tools.oneflow_gui.workflow_plan import build_case_graph_plan
from tools.oneflow_gui.workflow_process import run_process

DEFAULT_WORKFLOW = (
    REPOSITORY_ROOT / "grid/stl_to_gmsh/develop_test/stl_square_to_cfd.json"
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _load_workflow(path: Path) -> tuple[dict, list[dict]]:
    workflow = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(workflow, dict) or workflow.get("schema_version") != 1:
        raise ValueError("Unsupported development workflow definition.")
    steps = workflow.get("steps")
    if not isinstance(steps, list) or len(steps) != 1:
        raise ValueError("The workflow must contain exactly one case_graph step.")
    graph = steps[0]
    if not isinstance(graph, dict) or graph.get("type") != "case_graph":
        raise ValueError("The workflow step must be a case_graph.")
    cases = graph.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("The case_graph must contain a non-empty cases list.")
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("Every workflow case must be an object.")
        validation = case.get("validation")
        if validation is not None and not isinstance(validation, dict):
            raise ValueError(
                f"Case '{case.get('id', '<unknown>')}' validation must be an object."
            )
    return workflow, cases


def _ensure_outside_repository(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(REPOSITORY_ROOT)
    except ValueError:
        return resolved
    raise ValueError(f"Run data must stay outside the source repository: {resolved}")


def _case_file(case_directory: Path, value: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Validation artifact path must be a non-empty string.")
    relative = Path(value.strip())
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Validation artifact path must stay inside its Case.")
    return case_directory / relative


def _validate_step_outputs(case_directory: Path, specification: dict) -> dict:
    """Apply declarative, side-effect-free acceptance checks to one Case."""
    kind = specification.get("kind")
    if kind == "gmsh_msh4":
        mesh_path = _case_file(case_directory, specification.get("path"))
        mesh = inspect_msh4(mesh_path)
        expected_counts = specification.get("expected_counts")
        if not isinstance(expected_counts, dict):
            raise ValueError("gmsh_msh4 validation requires expected_counts.")
        observed_counts = {
            key: mesh[key]
            for key in (
                "elements",
                "nodes",
                "line_elements",
                "triangle_elements",
                "quadrilateral_elements",
            )
        }
        if observed_counts != expected_counts:
            raise ValueError(
                f"Gmsh mesh counts {observed_counts} do not match "
                f"expected {expected_counts}."
            )
        expected_element_types = {
            str(key): value
            for key, value in specification.get("expected_element_types", {}).items()
        }
        if mesh["element_types"] != expected_element_types:
            raise ValueError(
                f"Gmsh element types {mesh['element_types']} do not match "
                f"expected {expected_element_types}."
            )
        expected_groups = {
            (item["dimension"], item["name"])
            for item in specification.get("physical_groups", [])
        }
        actual_groups = {
            (item["dimension"], item["name"])
            for item in mesh["physical_groups"]
        }
        if actual_groups != expected_groups:
            raise ValueError(
                "Gmsh Physical Groups do not match the Case contract: "
                f"{sorted(actual_groups)} != {sorted(expected_groups)}."
            )
        return {
            "kind": kind,
            "path": str(mesh_path.relative_to(case_directory)),
            "mesh": observed_counts,
            "physical_groups": sorted(
                f"{dimension}:{name}" for dimension, name in actual_groups
            ),
        }

    if kind == "non_empty_files":
        paths = specification.get("paths")
        if not isinstance(paths, list) or not paths:
            raise ValueError("non_empty_files validation requires a non-empty paths list.")
        validated = []
        for relative in paths:
            path = _case_file(case_directory, relative)
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Required output is missing or empty: {path}")
            validated.append(
                {"path": str(path.relative_to(case_directory)), "size_bytes": path.stat().st_size}
            )
        return {"kind": kind, "files": validated}

    raise ValueError(f"Unsupported Case output validation kind: {kind!r}")


def run_pipeline(
    oneflow: Path,
    gmsh: Path,
    stage_root: Path,
    workflow_file: Path = DEFAULT_WORKFLOW,
) -> int:
    executable = oneflow.expanduser().resolve()
    gmsh_executable = gmsh.expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable}")
    if not gmsh_executable.is_file():
        raise FileNotFoundError(f"Gmsh executable not found: {gmsh_executable}")

    external_root = _ensure_outside_repository(stage_root)
    workflow_path = workflow_file.expanduser().resolve()
    workflow, cases = _load_workflow(workflow_path)
    plan = build_case_graph_plan(cases, REPOSITORY_ROOT)
    case_definitions = {case["id"]: case for case in cases}

    os.environ["PATH"] = (
        str(gmsh_executable.parent) + os.pathsep + os.environ.get("PATH", "")
    )

    run_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        + "-"
        + uuid4().hex[:8]
    )
    run_root = external_root / "pipelines" / run_id
    started_at = _utc_now()
    evidence = RunEnvironmentEvidence(
        run_id=run_id,
        executable=executable,
        run_directory=run_root,
        command_template=[str(executable), "r", "<case-directory>"],
        working_directory=REPOSITORY_ROOT,
        project_directory=REPOSITORY_ROOT / "grid/stl_to_gmsh",
        execution_mode="development_pipeline",
    )
    workflow_snapshot = run_root / "workflow.json"
    final_status = "failed"
    return_code = 1
    try:
        workflow_snapshot.write_text(
            workflow_path.read_text(encoding="utf-8"),
            encoding="utf-8",
            newline="\n",
        )

        def run_step(step):
            command = [str(executable), "r", str(step.directory)]
            output_lines: list[str] = []

            def on_output(line: str) -> None:
                output_lines.append(line)
                print(f"[{step.case_id}] {line}", flush=True)

            try:
                result = run_process(
                    command,
                    cwd=REPOSITORY_ROOT,
                    on_output=on_output,
                    is_cancelled=lambda: False,
                )
                status = result.status
                message = result.message
            except OSError as error:
                status = "failed"
                message = f"Unable to start OneFLOW: {error}"

            validation_path = None
            if status == "succeeded":
                validation_spec = case_definitions[step.case_id].get("validation")
                if validation_spec is not None:
                    try:
                        validation_report = _validate_step_outputs(
                            step.directory, validation_spec
                        )
                        validation_path = (
                            step.directory / "results/validation.json"
                        )
                        validation_path.parent.mkdir(parents=True, exist_ok=True)
                        validation_path.write_text(
                            json.dumps(validation_report, indent=2) + "\n",
                            encoding="utf-8",
                            newline="\n",
                        )
                    except (OSError, ValueError, KeyError, RuntimeError) as error:
                        status = "failed"
                        message = f"Artifact validation failed: {error}"

            results_directory = step.directory / "results"
            results_directory.mkdir(parents=True, exist_ok=True)
            log_path = results_directory / "development_case.log"
            log_path.write_text(
                "Command: "
                + json.dumps(command, ensure_ascii=True)
                + "\n"
                + "\n".join(output_lines)
                + ("\n" if output_lines else ""),
                encoding="utf-8",
                newline="\n",
            )
            print(f"[{step.case_id}] {status}: {message}", flush=True)
            return status, message

        graph_result = run_case_graph(plan, run_step, run_root=run_root)
        final_status = graph_result.status
        return_code = 0 if final_status == "succeeded" else 1

        step_records = []
        artifacts = []
        for result in graph_result.cases:
            staged_directory = run_root / "cases" / result.case_id
            plan_step = next(
                step for step in plan.steps if step.case_id == result.case_id
            )
            produced = []
            for output in plan_step.outputs:
                output_path = staged_directory / output.relative_path
                if output_path.is_file():
                    produced.append(str(output_path))
                    artifacts.append(str(output_path))
            validation_report = staged_directory / "results/validation.json"
            if validation_report.is_file():
                artifacts.append(str(validation_report))
            step_records.append(
                {
                    "case_id": result.case_id,
                    "status": result.status,
                    "message": result.message,
                    "directory": str(staged_directory),
                    "artifacts": produced,
                    "validation_report": (
                        str(validation_report) if validation_report.is_file() else None
                    ),
                    "log": str(staged_directory / "results/development_case.log"),
                }
            )

        report = {
            "schema_version": 1,
            "workflow_id": workflow["workflow_id"],
            "run_id": run_id,
            "status": final_status,
            "started_at": started_at,
            "finished_at": _utc_now(),
            "executable": str(executable),
            "gmsh_executable": str(gmsh_executable),
            "working_directory": str(REPOSITORY_ROOT),
            "workflow_file": str(workflow_path),
            "workflow_snapshot": str(workflow_snapshot),
            "run_directory": str(run_root),
            "steps": step_records,
            "artifacts": artifacts,
        }
        report_path = run_root / "pipeline-report.json"
        report_path.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        print(f"Pipeline report: {report_path}", flush=True)
        print(f"Pipeline status: {final_status}", flush=True)
        return return_code
    finally:
        evidence.finish(status=final_status, return_code=return_code)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--oneflow-executable", required=True, type=Path)
    parser.add_argument("--gmsh", required=True, type=Path)
    parser.add_argument(
        "--stage-root",
        type=Path,
        default=Path.home() / ".oneflow/mesh-layout-test/grid/stl_to_gmsh",
    )
    parser.add_argument("--workflow", type=Path, default=DEFAULT_WORKFLOW)
    args = parser.parse_args()
    try:
        return run_pipeline(
            args.oneflow_executable,
            args.gmsh,
            args.stage_root,
            args.workflow,
        )
    except (OSError, ValueError, KeyError, TypeError, RuntimeError, json.JSONDecodeError) as error:
        parser.exit(2, f"Development pipeline failed: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
