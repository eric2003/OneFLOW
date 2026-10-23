"""Execute a project or dependency workflow without a GUI process."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from .user_data_paths import UserDataPaths
from .workflow_graph_runner import run_case_graph
from .workflow_plan import build_case_graph_plan
from .workflow_process import run_process


def _log(message: str) -> None:
    print(message, flush=True)


def run_project(
    project_directory: Path,
    executable: Path,
    run_id: str,
    execution_mode: str,
    custom_root: Path | None,
    repository_root: Path,
) -> int:
    project_directory = project_directory.resolve()
    executable = executable.resolve()
    repository_root = repository_root.resolve()
    execution_directory = project_directory

    try:
        is_repository_content = (
            project_directory == repository_root
            or repository_root in project_directory.parents
        )
        if execution_mode == "source" and is_repository_content:
            raise ValueError(
                "Projects inside the OneFLOW source repository cannot run in place. "
                "Copy the project to a workspace first."
            )

        if execution_mode != "source":
            if execution_mode == "default_workspace":
                destination_root = UserDataPaths.discover().executions / run_id
            else:
                if custom_root is None:
                    raise ValueError("Choose a target directory for the project copy.")
                destination_root = custom_root
            destination_root = destination_root.resolve()
            if destination_root == repository_root or repository_root in destination_root.parents:
                raise ValueError("Project copies must be placed outside the OneFLOW source repository.")
            execution_directory = (destination_root / project_directory.name).resolve()
            if execution_directory == project_directory or project_directory in execution_directory.parents:
                raise ValueError("The project copy destination must be outside the source project.")
            destination_root.mkdir(parents=True, exist_ok=True)
            if execution_directory.exists():
                raise FileExistsError(f"Project copy already exists: {execution_directory}")
            shutil.copytree(
                project_directory,
                execution_directory,
                ignore=shutil.ignore_patterns(".git", "__pycache__"),
            )
            _log(f"Project copied to: {execution_directory}")

        workflow_file = execution_directory / "workflow.json"
        if workflow_file.is_file():
            workflow = json.loads(workflow_file.read_text(encoding="utf-8"))
            steps = workflow.get("steps")
            if not isinstance(steps, list):
                raise ValueError("Workflow definition must contain a steps list.")
            graph_steps = [
                step for step in steps
                if isinstance(step, dict) and step.get("type") == "case_graph"
            ]
            if len(graph_steps) != 1 or len(steps) != 1:
                raise ValueError("Project workflows must contain exactly one case_graph step.")
            plan = build_case_graph_plan(graph_steps[0].get("cases", []), execution_directory)
            run_root = execution_directory / ".oneflow-runs" / run_id
            _log(f"Workflow definition: {workflow_file}")
            _log(f"Workflow run directory: {run_root}")

            def run_step(step):
                command = [str(executable), "r", str(step.directory)]
                _log(f"Workflow case: {step.case_id}")
                _log("Command: " + json.dumps(command, ensure_ascii=True))
                result = run_process(command, cwd=executable.parent, on_output=_log)
                _log(f"Case {step.case_id}: {result.status} - {result.message}")
                return result.status, result.message

            graph_result = run_case_graph(
                plan,
                run_step,
                run_root=run_root,
                in_place=execution_mode == "source",
            )
            summary = "; ".join(
                f"{case.case_id}={case.status}" for case in graph_result.cases
            )
            _log(f"Workflow result: {summary}")
            return 0 if graph_result.status == "succeeded" else 1

        command = [str(executable), "r", str(execution_directory)]
        _log("Command: " + json.dumps(command, ensure_ascii=True))
        result = run_process(command, cwd=executable.parent, on_output=_log)
        _log(result.message)
        return result.return_code if result.return_code is not None else 1
    except (OSError, ValueError, shutil.Error, json.JSONDecodeError) as error:
        _log(f"Project run failed: {error}")
        return 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--execution-mode", required=True, choices=("source", "default_workspace", "custom_workspace"))
    parser.add_argument("--custom-root", type=Path)
    parser.add_argument("--repository-root", required=True, type=Path)
    args = parser.parse_args(argv)
    return run_project(
        args.project,
        args.executable,
        args.run_id,
        args.execution_mode,
        args.custom_root,
        args.repository_root,
    )


if __name__ == "__main__":
    raise SystemExit(main())
