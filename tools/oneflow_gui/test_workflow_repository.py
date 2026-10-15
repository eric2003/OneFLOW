"""Tests for workflow persistence that do not require Qt."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_repository import (
    WorkflowDefinition,
    WorkflowRepository,
    WorkflowRunRecord,
)


class WorkflowRepositoryTest(unittest.TestCase):
    def setUp(self) -> None:
        self._temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary_directory.cleanup)
        self.repository = WorkflowRepository(Path(self._temporary_directory.name))

    def test_definition_round_trips_and_can_be_updated(self) -> None:
        definition = WorkflowDefinition(
            workflow_id="mesh-check",
            name="Mesh check",
            description="Validate a mesh before solving.",
            steps=({"type": "validate_mesh", "strict": True},),
        )
        self.repository.save_workflow(definition)
        loaded = self.repository.load_workflow("mesh-check")
        self.assertEqual(loaded, definition)

        updated = WorkflowDefinition(
            workflow_id="mesh-check",
            name="Mesh check v2",
            steps=({"type": "validate_mesh"}, {"type": "generate_report"}),
        )
        self.repository.save_workflow(updated)
        self.assertEqual(self.repository.load_workflow("mesh-check"), updated)
        self.assertEqual(self.repository.list_workflows(), [updated])

    def test_unsupported_schema_version_is_rejected(self) -> None:
        path = Path(self._temporary_directory.name) / "definitions" / "future.json"
        path.parent.mkdir(parents=True)
        path.write_text(
            json.dumps({
                "schema_version": 99,
                "workflow_id": "future",
                "name": "Future",
                "steps": [],
            }),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "Unsupported workflow schema version"):
            self.repository.load_workflow("future")

    def test_workflow_deletion_preserves_run_history(self) -> None:
        self.repository.save_workflow(
            WorkflowDefinition("grid-run", "Grid run", steps=({"type": "run_oneflow"},))
        )
        self.repository.save_run(
            WorkflowRunRecord(
                run_id="run-001",
                workflow_id="grid-run",
                status="succeeded",
                started_at="2026-10-15T10:00:00Z",
                finished_at="2026-10-15T10:01:00Z",
                artifacts=("results/mesh-summary.json",),
            )
        )

        self.assertTrue(self.repository.delete_workflow("grid-run"))
        self.assertFalse(self.repository.delete_workflow("grid-run"))
        history = self.repository.list_runs("grid-run")
        self.assertEqual(len(history), 1)
        self.assertEqual(history[0].run_id, "run-001")
        self.assertEqual(history[0].artifacts, ("results/mesh-summary.json",))

    def test_run_history_is_filtered_and_sorted_newest_first(self) -> None:
        for run_id, workflow_id, started_at in (
            ("run-old", "workflow-a", "2026-10-14T10:00:00Z"),
            ("run-new", "workflow-a", "2026-10-15T10:00:00Z"),
            ("run-other", "workflow-b", "2026-10-16T10:00:00Z"),
        ):
            self.repository.save_run(
                WorkflowRunRecord(
                    run_id=run_id,
                    workflow_id=workflow_id,
                    status="succeeded",
                    started_at=started_at,
                )
            )
        self.assertEqual(
            [item.run_id for item in self.repository.list_runs("workflow-a")],
            ["run-new", "run-old"],
        )
        self.assertEqual(
            [item.run_id for item in self.repository.list_runs()],
            ["run-other", "run-new", "run-old"],
        )

    def test_invalid_identifier_and_step_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "workflow_id"):
            WorkflowDefinition("../escape", "Bad", steps=()).to_dict()
            self.repository.save_workflow(
                WorkflowDefinition("../escape", "Bad", steps=())
            )
        with self.assertRaisesRegex(ValueError, "non-empty type"):
            self.repository.save_workflow(
                WorkflowDefinition("bad-step", "Bad step", steps=({},))
            )


if __name__ == "__main__":
    unittest.main()
