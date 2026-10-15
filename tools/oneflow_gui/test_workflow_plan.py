"""Tests for dependency-aware execution plans and artifact paths."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_plan import build_case_graph_plan


class WorkflowPlanTest(unittest.TestCase):
    def test_grid_results_artifact_binds_solver_to_producer_case(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = [
                {
                    "id": "rae2822_roe_sa",
                    "directory": "test/cfd/rae2822_roe_sa",
                    "depends_on": ["rae2822grid"],
                    "inputs": [
                        {"artifact": "rae2822_grid", "setting": "gridFileName", "parameter_file": "script/cfd.txt"}
                    ],
                },
                {
                    "id": "rae2822grid",
                    "directory": "test/grid/rae2822convert",
                    "outputs": [
                        {"id": "rae2822_grid", "path": "results/mesh/final/rae2822_2d.ofl", "parameter_file": "script/grid.txt", "setting": "targetGridFileName"}
                    ],
                },
            ]
            plan = build_case_graph_plan(cases, root)

            self.assertEqual(
                [step.case_id for step in plan.steps],
                ["rae2822grid", "rae2822_roe_sa"],
            )
            producer = plan.steps[0]
            self.assertEqual(producer.outputs[0].relative_path, Path("results/mesh/final/rae2822_2d.ofl"))
            self.assertEqual(producer.outputs[0].parameter_file, Path("script/grid.txt"))
            self.assertEqual(producer.outputs[0].producer_setting, "targetGridFileName")
            solver = plan.steps[1]
            binding = solver.inputs[0]
            self.assertEqual(binding.consumer_setting, "gridFileName")
            self.assertEqual(
                binding.runtime_overrides,
                {
                    "sourceGridCaseDir": str((root / "test/grid/rae2822convert").resolve()),
                    "gridFileName": "results/rae2822_2d.ofl",
                },
            )
            self.assertEqual(
                binding.source_case_directory,
                (root / "test/grid/rae2822convert").resolve(),
            )
            self.assertEqual(
                binding.source_path,
                (root / "test/grid/rae2822convert/results/mesh/final/rae2822_2d.ofl").resolve(),
            )

    def test_grid_subdirectory_is_also_a_valid_artifact_location(self) -> None:
        plan = build_case_graph_plan(
            [
                {
                    "id": "mesh",
                    "directory": "grid/rae2822convert",
                    "outputs": [{"id": "grid", "path": "grid/rae2822_2d.ofl"}],
                },
                {
                    "id": "solver",
                    "directory": "cfd/rae2822",
                    "depends_on": ["mesh"],
                    "inputs": [{"artifact": "grid", "setting": "gridFileName", "parameter_file": "script/cfd.txt"}],
                },
            ],
            ".",
        )
        self.assertEqual(
            plan.steps[1].inputs[0].relative_path,
            Path("grid/rae2822_2d.ofl"),
        )

    def test_rejects_input_without_direct_dependency(self) -> None:
        with self.assertRaisesRegex(ValueError, "must depend directly"):
            build_case_graph_plan(
                [
                    {
                        "id": "mesh",
                        "directory": "mesh",
                        "outputs": [{"id": "grid", "path": "results/grid.ofl"}],
                    },
                    {
                        "id": "solver",
                        "directory": "solver",
                        "inputs": [{"artifact": "grid", "setting": "gridFileName", "parameter_file": "script/cfd.txt"}],
                    },
                ],
                ".",
            )

    def test_rejects_paths_escaping_producer_directory(self) -> None:
        with self.assertRaisesRegex(ValueError, "stay within its producer case"):
            build_case_graph_plan(
                [
                    {
                        "id": "mesh",
                        "directory": "mesh",
                        "outputs": [{"id": "grid", "path": "../grid.ofl"}],
                    }
                ],
                ".",
            )


if __name__ == "__main__":
    unittest.main()
