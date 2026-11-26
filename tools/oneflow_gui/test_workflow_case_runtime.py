"""Tests for isolated workflow case preparation."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.path_utils import paths_equal
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

    def test_reuses_existing_stage_and_preserves_open_output_directory(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mesh" / "script").mkdir(parents=True)
            (root / "mesh" / "script" / "grid.txt").write_text(
                'string targetGridFileName = "results/grid.msh";\n', encoding="utf-8"
            )
            (root / "mesh" / "results").mkdir()
            (root / "mesh" / "results" / "grid.msh").write_text("template", encoding="utf-8")
            plan = build_case_graph_plan([
                {"id": "mesh", "directory": "mesh", "outputs": [
                    {"id": "grid", "path": "results/grid.msh",
                     "parameter_file": "script/grid.txt", "setting": "targetGridFileName"},
                    {"id": "visual", "path": "results/visual"},
                ]},
            ], root)
            run_root = root / "validation"
            staged_case = run_root / "cases" / "mesh"
            open_visual = staged_case / "results" / "visual"
            open_visual.mkdir(parents=True)
            (open_visual / "view.dat").write_text("keep", encoding="utf-8")

            prepared = prepare_case_directory(plan.steps[0], run_root, {})

            self.assertTrue(paths_equal(prepared, staged_case))
            self.assertTrue((open_visual / "view.dat").is_file())
            self.assertFalse((staged_case / "results" / "grid.msh").exists())


    def test_stages_shared_parent_level_stl_cut_configuration(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case = root / "oblique_cut" / "cpp"
            (case / "script").mkdir(parents=True)
            (case / "input").mkdir()
            shared = root / "oblique_cut" / "script" / "cutting.json"
            shared.parent.mkdir(parents=True, exist_ok=True)
            shared.write_text('{"cutting_planes": {}}\n', encoding="utf-8")
            (case / "script" / "grid.txt").write_text(
                'string stlSurfaceCutConfigFile = "../script/cutting.json";\n'
                'string targetGridFileName = "results/cut.msh";\n',
                encoding="utf-8",
            )
            plan = build_case_graph_plan([
                {
                    "id": "oblique_cut",
                    "directory": "oblique_cut/cpp",
                    "outputs": [{
                        "id": "mesh",
                        "path": "results/cut.msh",
                        "parameter_file": "script/grid.txt",
                        "setting": "targetGridFileName",
                    }],
                },
            ], root)
            run_root = root / "validation"
            stale_config = run_root / "cases" / "script" / "cutting.json"
            stale_config.parent.mkdir(parents=True)
            stale_config.write_text('{"stale": true}\\n', encoding="utf-8")
            staged = prepare_case_directory(plan.steps[0], run_root, {})
            staged_shared = run_root / "cases" / "script" / "cutting.json"
            self.assertEqual(staged_shared.read_text(encoding="utf-8"), '{"cutting_planes": {}}\n')
            self.assertEqual(
                (staged / "script" / "grid.txt").read_text(encoding="utf-8").splitlines()[0],
                'string stlSurfaceCutConfigFile = "../script/cutting.json";',
            )
            self.assertEqual(shared.read_text(encoding="utf-8"), '{"cutting_planes": {}}\n')

    def test_stages_shared_cut_config_for_python_reference_case(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            case = root / "oblique_cut" / "python"
            (case / "input").mkdir(parents=True)
            shared = root / "oblique_cut" / "script" / "cutting.json"
            shared.parent.mkdir(parents=True)
            shared.write_text('{"shared": true}\n', encoding="utf-8")
            plan = build_case_graph_plan([
                {
                    "id": "oblique_python",
                    "directory": "oblique_cut/python",
                    "runtime_directory": "stl_to_gmsh",
                    "outputs": [{"id": "mesh", "path": "results/cut.msh"}],
                },
            ], root)
            run_root = root / "validation"
            staged = prepare_case_directory(plan.steps[0], run_root, {})
            self.assertTrue((run_root / "script" / "cutting.json").is_file())
            self.assertEqual(
                (run_root / "script" / "cutting.json").read_text(encoding="utf-8"),
                '{"shared": true}\n',
            )

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
            config.write_text('string targetGridFileName = "grid/old.ofl";\n', encoding="utf-8")
            stale_output = root / "mesh" / "results" / "generated" / "final" / "mesh.ofl"
            stale_output.parent.mkdir(parents=True)
            stale_output.write_bytes(b"stale")
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
            self.assertFalse((staged / "results" / "generated" / "final" / "mesh.ofl").exists())
            self.assertEqual(config.read_text(encoding="utf-8"), 'string targetGridFileName = "grid/old.ofl";\n')


if __name__ == "__main__":
    unittest.main()
