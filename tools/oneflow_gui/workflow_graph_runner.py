"""Execute a workflow plan in isolated case directories."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from pathlib import Path

from .workflow_case_runtime import prepare_case_directory
from .workflow_plan import WorkflowExecutionPlan, WorkflowPlanStep


@dataclass(frozen=True, slots=True)
class CaseExecutionResult:
    case_id: str
    status: str
    message: str


@dataclass(frozen=True, slots=True)
class WorkflowGraphRunResult:
    status: str
    cases: tuple[CaseExecutionResult, ...]


def run_case_graph(
    plan: WorkflowExecutionPlan,
    run_step: Callable[[WorkflowPlanStep], tuple[str, str]],
    *,
    run_root: Path | str,
    is_cancelled: Callable[[], bool] = lambda: False,
) -> WorkflowGraphRunResult:
    """Stage and run cases in dependency order.

    The callback receives a step whose directory and input producer directories
    point into this run's isolated workspace. It returns (status, message).
    """
    results: list[CaseExecutionResult] = []
    status_by_id: dict[str, str] = {}
    runtime_directories: dict[str, Path] = {}
    cancellation_requested = False

    for step in plan.steps:
        if cancellation_requested or is_cancelled():
            results.append(CaseExecutionResult(step.case_id, "cancelled", "Run cancelled."))
            status_by_id[step.case_id] = "cancelled"
            continue

        failed_dependencies = [
            dependency for dependency in step.depends_on
            if status_by_id.get(dependency) != "succeeded"
        ]
        if failed_dependencies:
            message = "Blocked by unsuccessful dependencies: " + ", ".join(failed_dependencies)
            results.append(CaseExecutionResult(step.case_id, "skipped", message))
            status_by_id[step.case_id] = "skipped"
            continue

        missing_inputs: list[str] = []
        for binding in step.inputs:
            producer_directory = runtime_directories.get(binding.producer_id)
            if producer_directory is None:
                missing_inputs.append(str(binding.source_path))
            else:
                artifact_path = producer_directory / binding.relative_path
                if not artifact_path.is_file():
                    missing_inputs.append(str(artifact_path))
        if missing_inputs:
            message = "Required input artifacts are missing: " + ", ".join(missing_inputs)
            results.append(CaseExecutionResult(step.case_id, "failed", message))
            status_by_id[step.case_id] = "failed"
            continue

        try:
            staged_directory = prepare_case_directory(
                step, run_root, runtime_directories
            )
        except (OSError, ValueError) as error:
            results.append(CaseExecutionResult(step.case_id, "failed", str(error)))
            status_by_id[step.case_id] = "failed"
            continue

        runtime_directories[step.case_id] = staged_directory
        staged_inputs = tuple(
            replace(
                binding,
                producer_directory=runtime_directories[binding.producer_id],
            )
            for binding in step.inputs
        )
        staged_step = replace(
            step,
            directory=staged_directory,
            inputs=staged_inputs,
        )

        status, message = run_step(staged_step)
        if status not in {"succeeded", "failed", "cancelled"}:
            raise ValueError(f"Unsupported case execution status: {status!r}")
        if status == "succeeded":
            missing_outputs = [
                str(staged_directory / output.relative_path)
                for output in step.outputs
                if not (staged_directory / output.relative_path).is_file()
            ]
            if missing_outputs:
                status = "failed"
                message = "Declared output artifacts were not produced: " + ", ".join(missing_outputs)
        results.append(CaseExecutionResult(step.case_id, status, message))
        status_by_id[step.case_id] = status
        if status == "cancelled":
            cancellation_requested = True

    if any(item.status == "failed" for item in results):
        overall_status = "failed"
    elif any(item.status == "cancelled" for item in results):
        overall_status = "cancelled"
    elif any(item.status == "skipped" for item in results):
        overall_status = "failed"
    else:
        overall_status = "succeeded"
    return WorkflowGraphRunResult(overall_status, tuple(results))
