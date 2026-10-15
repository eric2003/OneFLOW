"""Execute a workflow plan sequentially with dependency and artifact checks."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

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
    is_cancelled: Callable[[], bool] = lambda: False,
) -> WorkflowGraphRunResult:
    """Run planned cases in order and prevent consumers from using missing output.

    run_step returns (status, message). It is responsible for launching the
    case with the step's runtime input overrides, without editing source files.
    """
    results: list[CaseExecutionResult] = []
    status_by_id: dict[str, str] = {}

    for step in plan.steps:
        if is_cancelled():
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

        missing_inputs = [
            str(binding.source_path) for binding in step.inputs
            if not binding.source_path.is_file()
        ]
        if missing_inputs:
            message = "Required input artifacts are missing: " + ", ".join(missing_inputs)
            results.append(CaseExecutionResult(step.case_id, "failed", message))
            status_by_id[step.case_id] = "failed"
            continue

        status, message = run_step(step)
        if status not in {"succeeded", "failed", "cancelled"}:
            raise ValueError(f"Unsupported case execution status: {status!r}")
        results.append(CaseExecutionResult(step.case_id, status, message))
        status_by_id[step.case_id] = status

    if any(item.status == "failed" for item in results):
        overall_status = "failed"
    elif any(item.status == "cancelled" for item in results):
        overall_status = "cancelled"
    elif any(item.status == "skipped" for item in results):
        overall_status = "failed"
    else:
        overall_status = "succeeded"
    return WorkflowGraphRunResult(overall_status, tuple(results))
