"""Persistent registry for detached OneFLOW background tasks."""

from __future__ import annotations

import json
import os
import tempfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path


ACTIVE_STATUSES = frozenset({"starting", "running", "stopping"})
TERMINAL_STATUSES = frozenset({"succeeded", "failed", "stopped", "lost"})


def utc_now() -> str:
    """Return an ISO-8601 UTC timestamp."""
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True, slots=True)
class BackgroundTask:
    """Persisted metadata needed to rediscover and manage one background task."""

    task_id: str
    project_directory: str
    executable: str
    execution_mode: str
    host_pid: int | None
    process_pid: int | None
    status: str
    started_at: str
    updated_at: str
    log_file: str
    message: str = ""

    @classmethod
    def create(
        cls,
        *,
        task_id: str,
        project_directory: Path,
        executable: Path,
        execution_mode: str,
        log_file: Path,
    ) -> "BackgroundTask":
        now = utc_now()
        return cls(
            task_id=task_id,
            project_directory=str(project_directory.resolve()),
            executable=str(executable.resolve()),
            execution_mode=execution_mode,
            host_pid=None,
            process_pid=None,
            status="starting",
            started_at=now,
            updated_at=now,
            log_file=str(log_file.resolve()),
        )


class BackgroundTaskRegistry:
    """Store task records outside the source tree using atomic JSON replacement."""

    SCHEMA_VERSION = 1

    def __init__(self, root: Path) -> None:
        self.root = root
        self.records_directory = root / "background-tasks"
        self.records_directory.mkdir(parents=True, exist_ok=True)

    def _path(self, task_id: str) -> Path:
        if not task_id or Path(task_id).name != task_id or task_id in {".", ".."}:
            raise ValueError("Invalid background task ID.")
        return self.records_directory / f"{task_id}.json"

    def save(self, task: BackgroundTask) -> None:
        """Atomically persist a task record so restarts never read a partial JSON file."""
        path = self._path(task.task_id)
        payload = {
            "schema_version": self.SCHEMA_VERSION,
            "task": asdict(task),
        }
        fd, temporary_name = tempfile.mkstemp(
            prefix=f".{task.task_id}.", suffix=".tmp", dir=self.records_directory
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as stream:
                json.dump(payload, stream, ensure_ascii=True, indent=2)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary_name, path)
        finally:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass

    def load(self, task_id: str) -> BackgroundTask:
        """Load and validate a single persisted task."""
        path = self._path(task_id)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema_version") != self.SCHEMA_VERSION:
            raise ValueError(f"Unsupported background task schema in {path}.")
        task = BackgroundTask(**payload["task"])
        if task.task_id != task_id:
            raise ValueError(f"Background task ID does not match {path}.")
        return task

    def list_tasks(self) -> list[BackgroundTask]:
        """Return task records in stable start-time order."""
        tasks = []
        for path in self.records_directory.glob("*.json"):
            try:
                tasks.append(self.load(path.stem))
            except (OSError, ValueError, TypeError, KeyError, json.JSONDecodeError):
                continue
        return sorted(tasks, key=lambda task: (task.started_at, task.task_id))

    def list_active_tasks(self) -> list[BackgroundTask]:
        """Return only tasks that still require active management."""
        return [task for task in self.list_tasks() if task.status in ACTIVE_STATUSES]

    def remove(self, task_id: str) -> None:
        """Remove a record after its task has been deliberately forgotten."""
        try:
            self._path(task_id).unlink()
        except FileNotFoundError:
            pass
