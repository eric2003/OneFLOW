"""Tests translating portable workspaces to workflow execution drafts."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workspace_repository import SCHEMA_VERSION, WorkspaceDefinition
from tools.oneflow_gui.workspace_workflow_adapter import (
    build_workflow_definition_from_workspace,
    resolve_workspace_case_directory,
    store_workspace_case_path,
)


def make_workspace(mode: str, cases: list[dict]) -> WorkspaceDefinition:
    return WorkspaceDefinition.from_dict({
        "schema_version": SCHEMA_VERSION,
        "workspace_id": "sample-workspace",
        "name": "Sample workspace",
        "description": "Workspace description",
        "execution_mode": mode,
        "cases": cases,
    })


class WorkspaceWorkflowAdapterTest(unittest.TestCase):
    def test_relative_case_path_round_trips_against_workspace_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            workspace_path = Path(temporary_directory) / "nested" / "sample.ofw"
            resolved = resolve_workspace_case_directory("cases/mesh", workspace_path)
            stored = store_workspace_case_path("cases/mesh", workspace_path)
            loaded = resolve_workspace_case_directory(stored, workspace_path)

        self.assertEqual(stored, "cases/mesh")
        self.assertEqual(loaded, resolved)

    def test_external_case_path_remains_absolute(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            workspace_path = root / "workspace" / "sample.ofw"
            outside_case = root / "outside-case"
            stored = store_workspace_case_path(str(outside_case), workspace_path)

        self.assertEqual(stored, str(outside_case.resolve()))

    def test_single_case_becomes_one_directory_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            workspace_path = Path(temporary_directory) / "sample.ofw"
            workspace = make_workspace("single_case", [
                {"case_id": "one", "path": "cases/one"},
            ])
            workflow = build_workflow_definition_from_workspace(workspace, workspace_path)

        self.assertEqual(workflow.steps[0]["type"], "multi_directory")
        self.assertEqual(
            workflow.steps[0]["directories"],
            [str((workspace_path.parent / "cases" / "one").resolve())],
        )

    def test_sequential_workspace_preserves_directory_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            workspace_path = Path(temporary_directory) / "sample.ofw"
            workspace = make_workspace("sequential_cases", [
                {"case_id": "first", "path": "cases/first"},
                {"case_id": "second", "path": "cases/second"},
            ])
            workflow = build_workflow_definition_from_workspace(workspace, workspace_path)

        self.assertEqual(
            workflow.steps[0]["directories"],
            [
                str((workspace_path.parent / "cases" / "first").resolve()),
                str((workspace_path.parent / "cases" / "second").resolve()),
            ],
        )

    def test_independent_cases_use_existing_single_process_workflow(self) -> None:
        workspace = make_workspace("independent_cases", [
            {"case_id": "one", "path": "one"},
            {"case_id": "two", "path": "two"},
        ])
        workflow = build_workflow_definition_from_workspace(
            workspace, Path("workspace.ofw")
        )

        self.assertEqual(workflow.steps[0]["type"], "multi_directory")
        self.assertEqual(len(workflow.steps[0]["directories"]), 2)

    def test_dependency_workspace_preserves_graph_edges_and_base(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            workspace_path = Path(temporary_directory) / "sample.ofw"
            workspace = make_workspace("dependency_workflow", [
                {
                    "case_id": "producer",
                    "path": "cases/producer",
                    "outputs": [{
                        "id": "mesh",
                        "path": "results/mesh.ofl",
                        "parameter_file": "script/grid.txt",
                        "setting": "targetGridFileName",
                    }],
                },
                {
                    "case_id": "consumer",
                    "path": "cases/consumer",
                    "depends_on": ["producer"],
                    "inputs": [{
                        "artifact": "mesh",
                        "setting": "gridFileName",
                        "parameter_file": "script/cfd.txt",
                    }],
                },
            ])
            workflow = build_workflow_definition_from_workspace(workspace, workspace_path)

        step = workflow.steps[0]
        self.assertEqual(step["type"], "case_graph")
        self.assertEqual(step["base_directory"], str(workspace_path.parent.resolve()))
        self.assertEqual(step["cases"][1]["depends_on"], ["producer"])
        self.assertEqual(step["cases"][1]["directory"], "cases/consumer")
        self.assertEqual(step["cases"][0]["outputs"][0]["id"], "mesh")
        self.assertEqual(step["cases"][1]["inputs"][0]["artifact"], "mesh")

    def test_dependency_input_requires_declared_output(self) -> None:
        with self.assertRaisesRegex(ValueError, "unknown artifact"):
            make_workspace("dependency_workflow", [
                {"case_id": "producer", "path": "producer"},
                {
                    "case_id": "consumer",
                    "path": "consumer",
                    "depends_on": ["producer"],
                    "inputs": [{
                        "artifact": "missing",
                        "setting": "gridFileName",
                        "parameter_file": "script/cfd.txt",
                    }],
                },
            ])

    def test_artifact_producer_must_be_a_direct_dependency(self) -> None:
        with self.assertRaisesRegex(ValueError, "depend directly"):
            make_workspace("dependency_workflow", [
                {
                    "case_id": "producer",
                    "path": "producer",
                    "outputs": [{"id": "mesh", "path": "mesh.ofl"}],
                },
                {"case_id": "middle", "path": "middle", "depends_on": ["producer"]},
                {
                    "case_id": "consumer",
                    "path": "consumer",
                    "depends_on": ["middle"],
                    "inputs": [{
                        "artifact": "mesh",
                        "setting": "gridFileName",
                        "parameter_file": "script/cfd.txt",
                    }],
                },
            ])

    def test_artifact_setting_names_must_be_identifiers(self) -> None:
        with self.assertRaisesRegex(ValueError, "valid parameter setting name"):
            make_workspace("dependency_workflow", [
                {
                    "case_id": "producer",
                    "path": "producer",
                    "outputs": [{
                        "id": "mesh",
                        "path": "mesh.ofl",
                        "parameter_file": "grid.txt",
                        "setting": "invalid-name",
                    }],
                },
            ])

    def test_artifact_bindings_are_only_valid_for_dependency_workflows(self) -> None:
        with self.assertRaisesRegex(ValueError, "only in dependency workflows"):
            make_workspace("sequential_cases", [
                {
                    "case_id": "one",
                    "path": "one",
                    "outputs": [{"id": "mesh", "path": "mesh.ofl"}],
                },
            ])


if __name__ == "__main__":
    unittest.main()
