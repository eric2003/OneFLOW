#!/usr/bin/env python3
"""Run the staged STL -> Gmsh -> CGNS -> OneFLOW -> CFD development workflow."""

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
        raise ValueError("Unsupported STL-to-CFD workflow definition.")
    steps = workflow.get("steps")
    if not isinstance(steps, list) or len(steps) != 1:
        raise ValueError("The workflow must contain exactly one case_graph step.")
    graph = steps[0]
    if not isinstance(graph, dict) or graph.get("type") != "case_graph":
        raise ValueError("The workflow step must be a case_graph.")
    cases = graph.get("cases")
    if not isinstance(cases, list) or not cases:
        raise ValueError("The case_graph must contain a non-empty cases list.")
    if not isinstance(workflow.get("validation"), dict):
        raise ValueError("The workflow must declare artifact validation rules.")
    return workflow, cases


def _ensure_outside_repository(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    try:
        resolved.relative_to(REPOSITORY_ROOT)
    except ValueError:
        return resolved
    raise ValueError(f"Run data must stay outside the source repository: {resolved}")


def _validate_gmsh_mesh(mesh_path: Path, validation: dict) -> dict:
    mesh = inspect_msh4(mesh_path)
    expected = validation["gmsh_mesh"]
    observed = {
        key: mesh[key]
        for key in (
            "nodes",
            "line_elements",
            "triangle_elements",
            "quadrilateral_elements",
        )
    }
    if observed != expected:
        raise ValueError(
            f"Gmsh mesh counts {observed} do not match expected {expected}."
        )

    actual_groups = {
        (group["dimension"], group["name"])
        for group in mesh["physical_groups"]
    }
    expected_groups = {
        (group["dimension"], group["name"])
        for group in validation["physical_groups"]
    }
    if actual_groups != expected_groups:
        raise ValueError(
            "Gmsh Physical Groups do not match the pipeline contract: "
            f"{sorted(actual_groups)} != {sorted(expected_groups)}."
        )
    return {
        "mesh": observed,
        "physical_groups": sorted(
            f"{dimension}:{name}" for dimension, name in actual_groups
        ),
    }


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
        working_directory=executable.parent,
        project_directory=REPOSITORY_ROOT / "grid/stl_to_gmsh",
        execution_mode="development_pipeline",
    )

    final_status = "failed"
    return_code = 1
    try:
        def run_step(step):
            command = [str(executable), "r", str(step.directory)]
            output_lines: list[str] = []

            def on_output(line: str) -> None:
                output_lines.append(line)
                print(f"[{step.case_id}] {line}", flush=True)

            try:
                result = run_process(
                    command,
                    cwd=executable.parent,
                    on_output=on_output,
                    is_cancelled=lambda: False,
                )
                status = result.status
                message = result.message
            except OSError as error:
                status = "failed"
                message = f"Unable to start OneFLOW: {error}"

            if status == "succeeded":
                try:
                    if step.case_id == "stl_to_gmsh_square_2x2":
                        validation_report = _validate_gmsh_mesh(
                            step.directory / "results/square_2x2.msh",
                            workflow["validation"],
                        )
                        (step.directory / "results/mesh-validation.json").write_text(
                            json.dumps(validation_report, indent=2) + "\n",
                            encoding="utf-8",
                            newline="\n",
                        )
                    elif step.case_id == "cfd_solve_square_2x2":
                        output_paths = workflow["validation"]["cfd_non_empty_outputs"]
                        for relative in output_paths:
                            output = step.directory / relative
                            if not output.is_file() or output.stat().st_size == 0:
                                raise ValueError(
                                    f"CFD output is missing or empty: {output}"
                                )
                        (step.directory / "results/cfd-validation.json").write_text(
                            json.dumps(
                                {"non_empty_outputs": output_paths, "validated": True},
                                indent=2,
                            ) + "\n",
                            encoding="utf-8",
                            newline="\n",
                        )
                except (OSError, ValueError, KeyError) as error:
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
            step_records.append(
                {
                    "case_id": result.case_id,
                    "status": result.status,
                    "message": result.message,
                    "directory": str(staged_directory),
                    "artifacts": produced,
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
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
        parser.exit(2, f"STL-to-CFD pipeline failed: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
