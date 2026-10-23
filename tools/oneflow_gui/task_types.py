"""Shared task classifications and execution-location policies for the Workbench."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path


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


class TaskRunStatus(str, Enum):
    """Canonical persisted states shared by all Workbench task runners."""

    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"


_DISPLAY_STATUS_MAP = {
    "pending": TaskRunStatus.PENDING,
    "running": TaskRunStatus.RUNNING,
    "succeeded": TaskRunStatus.SUCCEEDED,
    "passed": TaskRunStatus.SUCCEEDED,
    "failed": TaskRunStatus.FAILED,
    "cancelled": TaskRunStatus.CANCELLED,
    "canceled": TaskRunStatus.CANCELLED,
}


def normalize_run_status(value: str | TaskRunStatus) -> TaskRunStatus:
    """Normalize persisted and GUI-facing status labels to one canonical state."""
    if isinstance(value, TaskRunStatus):
        return value
    if not isinstance(value, str):
        raise ValueError(f"Unsupported task run status: {value!r}.")
    status = _DISPLAY_STATUS_MAP.get(value.strip().lower())
    if status is None:
        raise ValueError(f"Unsupported task run status: {value!r}.")
    return status


@dataclass(frozen=True, slots=True)
class TaskRunResult:
    """Toolkit-independent outcome shared by Workbench task execution services."""

    task_type: TaskType
    status: TaskRunStatus
    return_code: int | None
    execution_directory: Path
    message: str

    def to_dict(self) -> dict[str, str | int | None]:
        """Return a JSON-compatible representation for run history."""
        return {
            "task_type": self.task_type.value,
            "status": self.status.value,
            "return_code": self.return_code,
            "execution_directory": str(self.execution_directory),
            "message": self.message,
        }

    def __post_init__(self) -> None:
        """Keep result values canonical at the boundary between runners and UI."""
        object.__setattr__(self, "task_type", TaskType(self.task_type))
        object.__setattr__(self, "status", normalize_run_status(self.status))
        object.__setattr__(
            self, "execution_directory", Path(self.execution_directory)
        )
        if self.return_code is not None and not isinstance(self.return_code, int):
            raise ValueError("Task result return_code must be an integer or null.")
        if not isinstance(self.message, str):
            raise ValueError("Task result message must be a string.")
