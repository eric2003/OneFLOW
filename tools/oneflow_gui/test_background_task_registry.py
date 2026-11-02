"""Tests for persistent background task registration."""

from __future__ import annotations

import json
import tempfile
from unittest.mock import patch
import unittest
from dataclasses import replace
from pathlib import Path

from tools.oneflow_gui.background_task_registry import (
    BackgroundTask,
    BackgroundTaskRegistry,
)


class BackgroundTaskRegistryTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.root = Path(self.temporary_directory.name)
        self.registry = BackgroundTaskRegistry(self.root)

    def make_task(self, task_id: str = "run-001") -> BackgroundTask:
        return BackgroundTask.create(
            task_id=task_id,
            project_directory=self.root / "project",
            executable=self.root / "OneFLOW.exe",
            execution_mode="source",
            log_file=self.root / "logs" / f"{task_id}.log",
        )

    def test_task_round_trips_across_registry_instances(self) -> None:
        task = self.make_task()
        self.registry.save(task)

        restarted_registry = BackgroundTaskRegistry(self.root)
        self.assertEqual(restarted_registry.load(task.task_id), task)
        self.assertEqual(restarted_registry.list_tasks(), [task])

    def test_save_replaces_record_without_leaving_temporary_files(self) -> None:
        task = self.make_task()
        self.registry.save(task)
        updated = replace(
            task,
            host_pid=1234,
            process_pid=5678,
            status="running",
            updated_at="2026-10-24T12:00:00+00:00",
        )
        self.registry.save(updated)

        self.assertEqual(self.registry.load(task.task_id), updated)
        self.assertEqual(
            [path.name for path in self.registry.records_directory.iterdir()],
            ["run-001.json"],
        )

    def test_save_retries_transient_windows_replace_denial(self) -> None:
        task = self.make_task()
        real_replace = __import__("os").replace
        attempts = 0

        def replace_with_transient_denial(source: str, destination: Path) -> None:
            nonlocal attempts
            attempts += 1
            if attempts < 3:
                raise PermissionError(5, "temporary sharing violation")
            real_replace(source, destination)

        with patch(
            "tools.oneflow_gui.background_task_registry.os.replace",
            side_effect=replace_with_transient_denial,
        ):
            self.registry.save(task)

        self.assertEqual(attempts, 3)
        self.assertEqual(self.registry.load(task.task_id), task)
        self.assertEqual(
            [path.name for path in self.registry.records_directory.iterdir()],
            [f"{task.task_id}.json"],
        )

    def test_save_preserves_existing_record_when_replace_remains_denied(self) -> None:
        original = self.make_task()
        self.registry.save(original)
        updated = replace(original, status="running")

        with patch(
            "tools.oneflow_gui.background_task_registry.os.replace",
            side_effect=PermissionError(5, "persistent access denial"),
        ):
            with self.assertRaises(PermissionError):
                self.registry.save(updated)

        self.assertEqual(self.registry.load(original.task_id), original)
        self.assertEqual(
            [path.name for path in self.registry.records_directory.iterdir()],
            [f"{original.task_id}.json"],
        )

    def test_invalid_task_id_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Invalid background task ID"):
            self.registry.load("../outside")

    def test_unknown_schema_version_is_rejected(self) -> None:
        task = self.make_task()
        self.registry.save(task)
        path = self.registry.records_directory / f"{task.task_id}.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["schema_version"] = 2
        path.write_text(json.dumps(payload), encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "Unsupported background task schema"):
            self.registry.load(task.task_id)

    def test_list_active_tasks_excludes_finished_records(self) -> None:
        starting = self.make_task("run-starting")
        running = replace(starting, task_id="run-running", status="running")
        finished = replace(starting, task_id="run-finished", status="succeeded")
        stopped = replace(starting, task_id="run-stopped", status="stopped")
        for task in (starting, running, finished, stopped):
            self.registry.save(task)

        self.assertEqual(
            {task.task_id for task in self.registry.list_active_tasks()},
            {"run-starting", "run-running"},
        )

    def test_remove_is_idempotent(self) -> None:
        self.registry.remove("run-001")
        self.registry.save(self.make_task())
        self.registry.remove("run-001")
        self.registry.remove("run-001")
        self.assertEqual(self.registry.list_tasks(), [])


if __name__ == "__main__":
    unittest.main()
