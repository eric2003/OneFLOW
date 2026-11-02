"""Tests for project file inventory, classification, and run changes."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.project_inventory import (
    compare_inventory,
    scan_project,
    scan_workflow_run_outputs,
)


class ProjectInventoryTest(unittest.TestCase):
    def test_scan_lists_directories_files_and_purpose_hints(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "script").mkdir()
            (root / "results").mkdir()
            (root / "workflow.json").write_text('{"steps": []}', encoding="utf-8")
            (root / "script" / "cfd.txt").write_text("string gridFileName = \"grid.ofl\";", encoding="utf-8")
            (root / "results" / "mesh.ofl").write_text("mesh", encoding="utf-8")

            entries = scan_project(root)

            self.assertEqual(entries["script"].kind, "directory")
            self.assertEqual(entries["workflow.json"].category, "Workflow definition")
            self.assertIn("workflow steps", entries["workflow.json"].purpose)
            self.assertEqual(entries["script/cfd.txt"].category, "Configuration")
            self.assertEqual(entries["results/mesh.ofl"].category, "Mesh data")

    def test_compare_inventory_reports_added_modified_deleted_and_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "same.txt").write_text("same", encoding="utf-8")
            (root / "change.txt").write_text("before", encoding="utf-8")
            (root / "delete.txt").write_text("delete", encoding="utf-8")
            before = scan_project(root)

            (root / "change.txt").write_text("after", encoding="utf-8")
            (root / "delete.txt").unlink()
            (root / "new.txt").write_text("new", encoding="utf-8")
            after = scan_project(root)

            changes = {change.entry.relative_path: change.status for change in compare_inventory(before, after)}

            self.assertEqual(changes["same.txt"], "unchanged")
            self.assertEqual(changes["change.txt"], "modified")
            self.assertEqual(changes["delete.txt"], "deleted")
            self.assertEqual(changes["new.txt"], "added")

    def test_scan_excludes_internal_run_staging_from_project_tree(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            staged = root / ".oneflow-runs" / "run-1" / "cases" / "mesh"
            staged.mkdir(parents=True)
            (staged / "mesh.ofl").write_text("mesh", encoding="utf-8")

            self.assertEqual(scan_project(root), {})

    def test_workflow_run_outputs_hide_unchanged_staged_inputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case = root / "mesh"
            (case / "grid").mkdir(parents=True)
            (case / "script").mkdir()
            (case / "grid" / "input.grd").write_text("input", encoding="utf-8")
            (case / "script" / "grid.txt").write_text("config", encoding="utf-8")
            (root / "workflow.json").write_text(
                '{"steps":[{"type":"case_graph","cases":[{"id":"mesh","directory":"mesh"}]}]}',
                encoding="utf-8",
            )
            run_case = root / ".oneflow-runs" / "run-1" / "cases" / "mesh"
            (run_case / "grid").mkdir(parents=True)
            (run_case / "script").mkdir()
            (run_case / "results").mkdir()
            (run_case / "grid" / "input.grd").write_text("input", encoding="utf-8")
            (run_case / "script" / "grid.txt").write_text("config", encoding="utf-8")
            (run_case / "results" / "mesh.ofl").write_text("generated", encoding="utf-8")

            outputs = scan_workflow_run_outputs(root, root / ".oneflow-runs" / "run-1")

            self.assertEqual(len(outputs), 1)
            self.assertIn("cases/mesh/results/mesh.ofl", outputs[0].relative_path)
            self.assertIn("Generated run file", outputs[0].purpose)


if __name__ == "__main__":
    unittest.main()
