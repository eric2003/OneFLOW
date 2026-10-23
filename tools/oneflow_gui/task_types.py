"""Shared task classifications and execution-location policies for the Workbench."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class TaskType(str, Enum):
    """User-visible categories of Workbench work."""

    CASE_DATABASE = "case_database"
    UNIT_TESTS = "unit_tests"
    CFD_GRID_REGRESSION = "cfd_grid_regression"


class ExecutionLocation(str, Enum):
    """Where a task is expected to execute."""

    SOURCE_BUILD = "source_build"
    USER_WORKSPACE = "user_workspace"


@dataclass(frozen=True, slots=True)
class TaskDefinition:
    """Stable task metadata shared by UI, execution, and history code."""

    task_type: TaskType
    title: str
    execution_location: ExecutionLocation
    supports_build_actions: bool = True


_TASK_DEFINITIONS = {
    TaskType.CASE_DATABASE: TaskDefinition(
        TaskType.CASE_DATABASE,
        "Case Database",
        ExecutionLocation.USER_WORKSPACE,
    ),
    TaskType.UNIT_TESTS: TaskDefinition(
        TaskType.UNIT_TESTS,
        "Unit Tests",
        ExecutionLocation.SOURCE_BUILD,
    ),
    TaskType.CFD_GRID_REGRESSION: TaskDefinition(
        TaskType.CFD_GRID_REGRESSION,
        "CFD/Grid Regression",
        ExecutionLocation.USER_WORKSPACE,
    ),
}


def task_definition(task_type: TaskType) -> TaskDefinition:
    """Return the shared definition for a task type."""
    return _TASK_DEFINITIONS[TaskType(task_type)]


def task_definitions() -> tuple[TaskDefinition, ...]:
    """Return task definitions in stable user-facing order."""
    return tuple(_TASK_DEFINITIONS.values())
