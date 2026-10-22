"""Tests for the portable OneFLOW workspace format."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workspace_repository import (
    WorkspaceDefinition,
    load_workspace,
    save_workspace,
)


def workspace_data(mode: str = "single_case", cases=None) -> dict:
    return {
        "schema_version": 1,
        "workspace_id": "sample-workspace",
        "name": "Sample workspace",
        "description": "",
        "execution_mode": mode,
        "cases": cases if cases is not None else [
            {"case_id": "case-a", "path": "cases/case-a", "depends_on": []}
        ],
    }


class WorkspaceRepositoryTest(unittest.TestCase):
    def test_round_trips_workspace_definition(self) -> None:
        workspace = WorkspaceDefinition.from_dict(workspace_data())
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = save_workspace(Path(temporary_directory) / "sample.ofw", workspace)
            loaded = load_workspace(path)
        self.assertEqual(loaded, workspace)

    def test_sequential_cases_preserve_declared_order(self) -> None:
        data = workspace_data("sequential_cases", [
            {"case_id": "first", "path": "first"},
            {"case_id": "second", "path": "second"},
        ])
        workspace = WorkspaceDefinition.from_dict(data)
        self.assertEqual([case.case_id for case in workspace.cases], ["first", "second"])

    def test_independent_cases_reject_dependencies(self) -> None:
        data = workspace_data("independent_cases", [
            {"case_id": "first", "path": "first"},
            {"case_id": "second", "path": "second", "depends_on": ["first"]},
        ])
        with self.assertRaisesRegex(ValueError, "cannot declare dependencies"):
            WorkspaceDefinition.from_dict(data)

    def test_dependency_workflow_rejects_cycles(self) -> None:
        data = workspace_data("dependency_workflow", [
            {"case_id": "first", "path": "first", "depends_on": ["second"]},
            {"case_id": "second", "path": "second", "depends_on": ["first"]},
        ])
        with self.assertRaisesRegex(ValueError, "contains a cycle"):
            WorkspaceDefinition.from_dict(data)

    def test_dependency_workflow_rejects_unknown_case(self) -> None:
        data = workspace_data("dependency_workflow", [
            {"case_id": "first", "path": "first", "depends_on": ["missing"]},
        ])
        with self.assertRaisesRegex(ValueError, "unknown case"):
            WorkspaceDefinition.from_dict(data)

    def test_single_case_requires_exactly_one_case(self) -> None:
        data = workspace_data("single_case", [
            {"case_id": "first", "path": "first"},
            {"case_id": "second", "path": "second"},
        ])
        with self.assertRaisesRegex(ValueError, "exactly one case"):
            WorkspaceDefinition.from_dict(data)

    def test_rejects_non_workspace_extension(self) -> None:
        workspace = WorkspaceDefinition.from_dict(workspace_data())
        with tempfile.TemporaryDirectory() as temporary_directory:
            with self.assertRaisesRegex(ValueError, "extension"):
                save_workspace(Path(temporary_directory) / "sample.json", workspace)


if __name__ == "__main__":
    unittest.main()
