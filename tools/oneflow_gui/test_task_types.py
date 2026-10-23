"""Tests for the shared Workbench task classifications."""

from __future__ import annotations

import unittest
from pathlib import Path

from tools.oneflow_gui.task_types import (
    ExecutionLocation,
    TaskRunResult,
    TaskRunStatus,
    TaskType,
    aggregate_run_status,
    normalize_run_status,
    task_definition,
    task_definitions,
)


class TaskTypesTest(unittest.TestCase):
    def test_task_types_keep_their_execution_locations(self) -> None:
        self.assertEqual(
            task_definition(TaskType.CASE_DATABASE).execution_location,
            ExecutionLocation.USER_WORKSPACE,
        )
        self.assertEqual(
            task_definition(TaskType.UNIT_TESTS).execution_location,
            ExecutionLocation.SOURCE_BUILD,
        )
        self.assertEqual(
            task_definition(TaskType.CFD_GRID_REGRESSION).execution_location,
            ExecutionLocation.USER_WORKSPACE,
        )

    def test_definitions_have_unique_task_types_and_titles(self) -> None:
        definitions = task_definitions()
        self.assertEqual(
            len({definition.task_type for definition in definitions}),
            len(definitions),
        )
        self.assertEqual(
            len({definition.title for definition in definitions}),
            len(definitions),
        )
        self.assertEqual(
            {definition.task_type for definition in definitions},
            set(TaskType),
        )

    def test_run_statuses_normalize_to_canonical_values(self) -> None:
        self.assertEqual(
            normalize_run_status("Pending"), TaskRunStatus.PENDING
        )
        self.assertEqual(
            normalize_run_status("Passed"), TaskRunStatus.SUCCEEDED
        )
        self.assertEqual(
            normalize_run_status("Succeeded"), TaskRunStatus.SUCCEEDED
        )
        self.assertEqual(
            normalize_run_status("Cancelled"), TaskRunStatus.CANCELLED
        )
        self.assertEqual(
            {status.value for status in TaskRunStatus},
            {"pending", "running", "succeeded", "failed", "cancelled"},
        )

    def test_queue_status_aggregates_success_failure_and_cancellation(self) -> None:
        self.assertEqual(
            aggregate_run_status(["succeeded", "succeeded"]),
            TaskRunStatus.SUCCEEDED,
        )
        self.assertEqual(
            aggregate_run_status(["succeeded", "cancelled"]),
            TaskRunStatus.CANCELLED,
        )
        self.assertEqual(
            aggregate_run_status(["cancelled", "failed", "succeeded"]),
            TaskRunStatus.FAILED,
        )
        self.assertEqual(
            aggregate_run_status(["failed", "cancelled"]),
            TaskRunStatus.FAILED,
        )

    def test_empty_queue_status_preserves_success_default(self) -> None:
        self.assertEqual(aggregate_run_status([]), TaskRunStatus.SUCCEEDED)

    def test_unknown_run_status_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unsupported task run status"):
            normalize_run_status("skipped")

    def test_task_run_result_normalizes_shared_values(self) -> None:
        result = TaskRunResult(
            task_type="unit_tests",
            status="Passed",
            return_code=0,
            execution_directory="build",
            message="All tests passed.",
        )
        self.assertEqual(result.task_type, TaskType.UNIT_TESTS)
        self.assertEqual(result.status, TaskRunStatus.SUCCEEDED)
        self.assertEqual(result.return_code, 0)
        self.assertEqual(result.execution_directory, Path("build"))
        self.assertEqual(result.message, "All tests passed.")

    def test_task_run_result_serializes_for_history(self) -> None:
        result = TaskRunResult(
            task_type=TaskType.CASE_DATABASE,
            status=TaskRunStatus.FAILED,
            return_code=2,
            execution_directory=Path("workspace/case-a"),
            message="case-a: Process failed.",
        )
        self.assertEqual(
            result.to_dict(),
            {
                "task_type": "case_database",
                "status": "failed",
                "return_code": 2,
                "execution_directory": "workspace/case-a",
                "message": "case-a: Process failed.",
            },
        )

    def test_task_run_result_allows_missing_return_code(self) -> None:
        result = TaskRunResult(
            task_type=TaskType.CFD_GRID_REGRESSION,
            status=TaskRunStatus.CANCELLED,
            return_code=None,
            execution_directory=Path("workspace"),
            message="Stopped by user.",
        )
        self.assertIsNone(result.return_code)
        self.assertEqual(result.status, TaskRunStatus.CANCELLED)

    def test_task_run_result_rejects_invalid_return_code(self) -> None:
        with self.assertRaisesRegex(ValueError, "return_code"):
            TaskRunResult(
                task_type=TaskType.UNIT_TESTS,
                status=TaskRunStatus.FAILED,
                return_code="1",
                execution_directory=Path("build"),
                message="Failed.",
            )

    def test_string_task_type_resolves_to_definition(self) -> None:
        definition = task_definition(TaskType("unit_tests"))
        self.assertEqual(definition.title, "Unit Tests")
        self.assertTrue(definition.supports_build_actions)


if __name__ == "__main__":
    unittest.main()
