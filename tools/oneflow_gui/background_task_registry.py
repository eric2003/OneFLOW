"""Persistent registry for detached OneFLOW background tasks."""

from __future__ import annotations

import json
import os
import tempfile
import time
from contextlib import contextmanager
from typing import Callable, Iterator
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
    REPLACE_RETRY_DELAYS = (0.01, 0.02, 0.04, 0.08, 0.16)

    def __init__(self, root: Path) -> None:
        self.root = root
        self.records_directory = root / "background-tasks"
        self.records_directory.mkdir(parents=True, exist_ok=True)

    def _path(self, task_id: str) -> Path:
        if not task_id or Path(task_id).name != task_id or task_id in {".", ".."}:
            raise ValueError("Invalid background task ID.")
        return self.records_directory / f"{task_id}.json"

    @contextmanager
    def _exclusive_update(self) -> Iterator[None]:
        """Serialize registry read-modify-write operations across processes."""
        lock_path = self.root / ".background-task-registry.lock"
        with lock_path.open("a+b") as stream:
            stream.seek(0, os.SEEK_END)
            if stream.tell() == 0:
                stream.write(b"\0")
                stream.flush()
            stream.seek(0)
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(stream.fileno(), msvcrt.LK_LOCK, 1)
            else:
                import fcntl

                fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                if os.name == "nt":
                    msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
                else:
                    fcntl.flock(stream.fileno(), fcntl.LOCK_UN)

    def _save_unlocked(self, task: BackgroundTask) -> None:
        """Atomically persist a task record while the caller owns the registry lock."""
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
            for attempt, delay in enumerate((0.0, *self.REPLACE_RETRY_DELAYS)):
                if delay:
                    time.sleep(delay)
                try:
                    os.replace(temporary_name, path)
                    break
                except PermissionError:
                    if attempt == len(self.REPLACE_RETRY_DELAYS):
                        raise
        finally:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass

    def save(self, task: BackgroundTask) -> None:
        """Serialize and atomically persist one task record."""
        with self._exclusive_update():
            self._save_unlocked(task)

    def _load_unlocked(self, task_id: str) -> BackgroundTask:
        """Load and validate a task while the caller owns the registry lock."""
        path = self._path(task_id)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema_version") != self.SCHEMA_VERSION:
            raise ValueError(f"Unsupported background task schema in {path}.")
        task = BackgroundTask(**payload["task"])
        if task.task_id != task_id:
            raise ValueError(f"Background task ID does not match {path}.")
        return task

    def load(self, task_id: str) -> BackgroundTask:
        """Load a complete record; atomic replacement protects readers."""
        return self._load_unlocked(task_id)

    def update(
        self,
        task_id: str,
        transform: Callable[[BackgroundTask], BackgroundTask],
    ) -> BackgroundTask:
        """Apply a read-modify-write transition without losing concurrent updates."""
        with self._exclusive_update():
            current = self._load_unlocked(task_id)
            updated = transform(current)
            if updated.task_id != task_id:
                raise ValueError("A background task update cannot change its ID.")
            self._save_unlocked(updated)
            return updated

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
        """Remove a record without racing a concurrent state transition."""
        with self._exclusive_update():
            try:
                self._path(task_id).unlink()
            except FileNotFoundError:
                pass
