"""Tests for isolated workflow case preparation."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_case_runtime import prepare_case_directory
from tools.oneflow_gui.workflow_plan import build_case_graph_plan


class WorkflowCaseRuntimeTest(unittest.TestCase):
    def test_stages_consumer_and_preserves_source_config(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mesh").mkdir()
            (root / "solver" / "script").mkdir(parents=True)
            original = root / "solver" / "script" / "cfd.txt"
            original.write_text('string gridFileName = "grid/old.ofl";\nint maxSteps = 50;\n', encoding="utf-8")
            plan = build_case_graph_plan([
                {"id": "mesh", "directory": "mesh", "outputs": [{"id": "grid", "path": "results/grid.ofl"}]},
                {"id": "solver", "directory": "solver", "depends_on": ["mesh"],
                 "inputs": [{"artifact": "grid", "setting": "gridFileName", "parameter_file": "script/cfd.txt"}]},
            ], root)
            producer_run = root / "run" / "cases" / "mesh"
            (producer_run / "results").mkdir(parents=True)
            (producer_run / "results" / "grid.ofl").write_bytes(b"grid")
            step = next(item for item in plan.steps if item.case_id == "solver")
            staged = prepare_case_directory(step, root / "run", {"mesh": producer_run})
            config = (staged / "script" / "cfd.txt").read_text(encoding="utf-8")
            escaped_producer_path = str(producer_run.resolve()).replace("\\", "\\\\")
            self.assertIn(f'string sourceGridCaseDir = "{escaped_producer_path}";', config)
            self.assertIn('string gridFileName = "results/grid.ofl";', config)
            self.assertEqual(original.read_text(encoding="utf-8"), 'string gridFileName = "grid/old.ofl";\nint maxSteps = 50;\n')

    def test_rejects_parameter_file_traversal(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mesh").mkdir()
            (root / "solver").mkdir()
            with self.assertRaisesRegex(ValueError, "stay within its producer case"):
                build_case_graph_plan([
                    {"id": "mesh", "directory": "mesh", "outputs": [{"id": "grid", "path": "grid.ofl"}]},
                    {"id": "solver", "directory": "solver", "depends_on": ["mesh"],
                     "inputs": [{"artifact": "grid", "setting": "gridFileName", "parameter_file": "../outside.txt"}]},
                ], root)


    def test_configures_producer_output_in_arbitrary_nested_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mesh" / "script").mkdir(parents=True)
            config = root / "mesh" / "script" / "grid.txt"
            config.write_text('string targetGridFileName = "grid/old.ofl";\\n', encoding="utf-8")
            plan = build_case_graph_plan([
                {
                    "id": "mesh",
                    "directory": "mesh",
                    "outputs": [{
                        "id": "grid",
                        "path": "results/generated/final/mesh.ofl",
                        "parameter_file": "script/grid.txt",
                        "setting": "targetGridFileName",
                    }],
                },
            ], root)
            staged = prepare_case_directory(plan.steps[0], root / "run", {})
            staged_config = (staged / "script" / "grid.txt").read_text(encoding="utf-8")
            self.assertIn('string targetGridFileName = "results/generated/final/mesh.ofl";', staged_config)
            self.assertTrue((staged / "results" / "generated" / "final").is_dir())
            self.assertEqual(config.read_text(encoding="utf-8"), 'string targetGridFileName = "grid/old.ofl";\\n')


if __name__ == "__main__":
    unittest.main()
