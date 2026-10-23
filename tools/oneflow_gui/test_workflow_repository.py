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

    def test_case_graph_round_trips_with_explicit_dependencies(self) -> None:
        definition = WorkflowDefinition(
            workflow_id="mesh-to-solver",
            name="Mesh to solver",
            steps=(
                {
                    "type": "case_graph",
                    "cases": [
                        {"id": "mesh", "directory": "grid/generated"},
                        {
                            "id": "solve",
                            "directory": "test/cfd/solver",
                            "depends_on": ["mesh"],
                        },
                    ],
                },
            ),
        )
        self.repository.save_workflow(definition)
        self.assertEqual(self.repository.load_workflow("mesh-to-solver"), definition)

    def test_invalid_case_graph_is_rejected_when_saving(self) -> None:
        definition = WorkflowDefinition(
            workflow_id="cyclic-workflow",
            name="Cyclic workflow",
            steps=(
                {
                    "type": "case_graph",
                    "cases": [
                        {"id": "mesh", "directory": "grid", "depends_on": ["solve"]},
                        {"id": "solve", "directory": "solver", "depends_on": ["mesh"]},
                    ],
                },
            ),
        )
        with self.assertRaisesRegex(ValueError, "contain a cycle"):
            self.repository.save_workflow(definition)

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

    def test_run_notes_and_record_deletion_do_not_touch_artifacts(self) -> None:
        self.repository.save_run(
            WorkflowRunRecord(
                run_id="run-note",
                workflow_id="removed-workflow",
                status="succeeded",
                started_at="2026-10-15T10:00:00Z",
                artifacts=("executions/run-note/results/mesh.ofl",),
            )
        )
        note_path = Path(self._temporary_directory.name) / "executions" / "run-note" / "results" / "mesh.ofl"
        note_path.parent.mkdir(parents=True)
        note_path.write_text("artifact", encoding="utf-8")

        updated = self.repository.update_run_note("run-note", "Reviewed mesh output.")
        self.assertEqual(updated.note, "Reviewed mesh output.")
        self.assertEqual(self.repository.list_runs()[0].note, "Reviewed mesh output.")
        self.assertTrue(self.repository.delete_run("run-note"))
        self.assertFalse(self.repository.delete_run("run-note"))
        self.assertTrue(note_path.is_file())
        self.assertEqual(note_path.read_text(encoding="utf-8"), "artifact")

    def test_history_cleanup_filters_by_status_workflow_and_age(self) -> None:
        records = (
            ("old-failed", "workflow-a", "failed", "2026-10-10T10:00:00Z"),
            ("old-passed", "workflow-a", "succeeded", "2026-10-11T10:00:00Z"),
            ("new-failed", "workflow-a", "failed", "2026-10-16T10:00:00Z"),
            ("other-failed", "workflow-b", "failed", "2026-10-09T10:00:00Z"),
        )
        for run_id, workflow_id, status, started_at in records:
            self.repository.save_run(
                WorkflowRunRecord(
                    run_id=run_id,
                    workflow_id=workflow_id,
                    status=status,
                    started_at=started_at,
                )
            )

        self.assertEqual(
            [record.run_id for record in self.repository.list_runs(
                "workflow-a", status="failed", started_before="2026-10-15T00:00:00Z"
            )],
            ["old-failed"],
        )
        removed = self.repository.clear_runs(
            workflow_id="workflow-a",
            statuses={"failed"},
            started_before="2026-10-15T00:00:00Z",
        )
        self.assertEqual(removed, 1)
        self.assertEqual(
            {record.run_id for record in self.repository.list_runs()},
            {"old-passed", "new-failed", "other-failed"},
        )

    def test_task_type_round_trips_and_legacy_records_default_to_workflow(self) -> None:
        record = WorkflowRunRecord(
            run_id="run-case-db",
            workflow_id="case-workbench",
            status="running",
            started_at="2026-10-23T10:00:00Z",
            task_type="case_database",
        )
        self.repository.save_run(record)
        loaded = self.repository.list_runs()[0]
        self.assertEqual(loaded.task_type, "case_database")

        legacy = WorkflowRunRecord.from_dict({
            "run_id": "legacy-run",
            "workflow_id": "old-workflow",
            "status": "succeeded",
            "started_at": "2026-10-22T10:00:00Z",
        })
        self.assertEqual(legacy.task_type, "workflow")

    def test_empty_task_type_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "task_type"):
            WorkflowRunRecord.from_dict({
                "run_id": "bad-task",
                "workflow_id": "workflow",
                "status": "running",
                "started_at": "2026-10-23T10:00:00Z",
                "task_type": "",
            })

    def test_run_log_path_round_trips(self) -> None:
        record = WorkflowRunRecord(
            run_id="run-log",
            workflow_id="workflow-log",
            status="succeeded",
            started_at="2026-10-17T10:00:00Z",
            log_path="/tmp/oneflow-run.log",
        )
        self.repository.save_run(record)
        self.assertEqual(self.repository.list_runs()[0].log_path, "/tmp/oneflow-run.log")

    def test_legacy_display_status_is_normalized_when_loaded(self) -> None:
        record = WorkflowRunRecord.from_dict({
            "run_id": "legacy-passed",
            "workflow_id": "old-task",
            "status": "Passed",
            "started_at": "2026-10-23T10:00:00Z",
        })
        self.assertEqual(record.status, "succeeded")

    def test_history_filters_reject_unknown_status(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unsupported workflow run status"):
            self.repository.list_runs(status="unknown")
        with self.assertRaisesRegex(ValueError, "Unsupported workflow run status"):
            self.repository.clear_runs(statuses={"unknown"})

    def test_invalid_identifier_and_step_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "workflow_id"):
            self.repository.save_workflow(
                WorkflowDefinition("../escape", "Bad", steps=())
            )
        with self.assertRaisesRegex(ValueError, "non-empty type"):
            self.repository.save_workflow(
                WorkflowDefinition("bad-step", "Bad step", steps=({},))
            )
        with self.assertRaisesRegex(ValueError, "run status"):
            self.repository.save_run(
                WorkflowRunRecord(
                    run_id="bad-status",
                    workflow_id="workflow",
                    status="unknown",
                    started_at="2026-10-15T10:00:00Z",
                )
            )


if __name__ == "__main__":
    unittest.main()
