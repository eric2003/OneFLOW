"""Small, immutable data models shared by the desktop application."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class RunState(str, Enum):
    """Lifecycle states displayed for each case in the run queue."""

    PENDING = "Pending"
    RUNNING = "Running"
    SUCCEEDED = "Succeeded"
    FAILED = "Failed"
    CANCELLED = "Cancelled"


@dataclass(frozen=True, slots=True)
class CaseDefinition:
    """Resolved catalog entry; paths are absolute and validated at load time."""

    case_id: str
    title: str
    case_directory: Path
    runtime_config: Path
    profile: Path | None
    layout: Path | None


@dataclass(frozen=True, slots=True)
class CaseRunResult:
    """Outcome of one process invocation, independent of the GUI toolkit."""

    case_id: str
    state: RunState
    return_code: int | None
    workspace: Path
    message: str
