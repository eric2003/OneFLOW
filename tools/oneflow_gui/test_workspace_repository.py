"""Tests for the portable OneFLOW workspace format."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workspace_repository import (
    SCHEMA_VERSION,
    WorkspaceDefinition,
    load_workspace,
    save_workspace,
)


def workspace_data(mode: str = "single_case", cases=None) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
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

    def test_version_one_workspace_migrates_without_losing_cases(self) -> None:
        data = workspace_data()
        data["schema_version"] = 1
        workspace = WorkspaceDefinition.from_dict(data)
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = save_workspace(
                Path(temporary_directory) / "migrated.ofw", workspace
            )
            persisted = json.loads(path.read_text(encoding="utf-8"))

        self.assertEqual(workspace.schema_version, SCHEMA_VERSION)
        self.assertEqual(persisted["schema_version"], SCHEMA_VERSION)
        self.assertEqual(workspace.cases[0].inputs, ())
        self.assertEqual(workspace.cases[0].outputs, ())

    def test_unknown_workspace_schema_version_is_rejected(self) -> None:
        data = workspace_data()
        data["schema_version"] = 99

        with self.assertRaisesRegex(ValueError, "Unsupported workspace schema version"):
            WorkspaceDefinition.from_dict(data)

    def test_version_one_cannot_contain_version_two_artifact_fields(self) -> None:
        data = workspace_data("dependency_workflow", [
            {"case_id": "grid", "path": "grid", "outputs": [{"id": "mesh", "path": "mesh.ofl"}]},
        ])
        data["schema_version"] = 1

        with self.assertRaisesRegex(ValueError, "require workspace schema version 2"):
            WorkspaceDefinition.from_dict(data)

    def test_dependency_artifact_bindings_round_trip(self) -> None:
        data = workspace_data("dependency_workflow", [
            {
                "case_id": "producer",
                "path": "producer",
                "outputs": [{
                    "id": "mesh",
                    "path": "results/mesh.ofl",
                    "parameter_file": "script/grid.txt",
                    "setting": "targetGridFileName",
                }],
            },
            {
                "case_id": "consumer",
                "path": "consumer",
                "depends_on": ["producer"],
                "inputs": [{
                    "artifact": "mesh",
                    "setting": "gridFileName",
                    "parameter_file": "script/cfd.txt",
                }],
            },
        ])
        workspace = WorkspaceDefinition.from_dict(data)
        with tempfile.TemporaryDirectory() as temporary_directory:
            path = save_workspace(Path(temporary_directory) / "artifacts.ofw", workspace)
            loaded = load_workspace(path)

        self.assertEqual(loaded, workspace)
        self.assertEqual(loaded.cases[0].outputs[0].artifact_id, "mesh")
        self.assertEqual(loaded.cases[1].inputs[0].artifact, "mesh")

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
