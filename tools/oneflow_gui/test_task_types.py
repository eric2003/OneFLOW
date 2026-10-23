"""Tests for the shared Workbench task classifications."""

from __future__ import annotations

import unittest

from tools.oneflow_gui.task_types import (
    ExecutionLocation,
    TaskRunStatus,
    TaskType,
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

    def test_unknown_run_status_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unsupported task run status"):
            normalize_run_status("skipped")

    def test_string_task_type_resolves_to_definition(self) -> None:
        definition = task_definition(TaskType("unit_tests"))
        self.assertEqual(definition.title, "Unit Tests")
        self.assertTrue(definition.supports_build_actions)


if __name__ == "__main__":
    unittest.main()
