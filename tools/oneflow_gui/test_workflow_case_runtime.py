"""Tests for isolated workflow case preparation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_plan import build_case_graph_plan
from tools.oneflow_gui.workflow_case_runtime import prepare_case_directory


class WorkflowCaseRuntimeTest(unittest.TestCase):
    def test_stages_consumer_and_binds_grid_from_staged_producer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            producer = root / "mesh"
            producer.mkdir()
            consumer = root / "solver"
            (consumer / "script").mkdir(parents=True)
            original_config = consumer / "script" / "cfd.txt"
            original_config.write_text(
                'string gridFileName = "grid/old.ofl";\n'
                'int maxSteps = 50;\n',
                encoding="utf-8",
            )
            plan = build_case_graph_plan(
                [
                    {
                        "id": "mesh",
                        "directory": "mesh",
                        "outputs": [{"id": "grid", "path": "results/grid.ofl"}],
                    },
                    {
                        "id": "solver",
                        "directory": "solver",
                        "depends_on": ["mesh"],
                        "inputs": [
                            {
                                "artifact": "grid",
                                "setting": "gridFileName",
                                "parameter_file": "script/cfd.txt",
                            }
                        ],
                    },
                ],
                root,
            )
            producer_run = root / "run" / "cases" / "mesh"
            (producer_run / "results").mkdir(parents=True)
            (producer_run / "results" / "grid.ofl").write_bytes(b"mesh")
            solver_step = next(step for step in plan.steps if step.case_id == "solver")

            staged = prepare_case_directory(
                solver_step,
                root / "run",
                {"mesh": producer_run},
            )

            config = (staged / "script" / "cfd.txt").read_text(encoding="utf-8")
            self.assertIn(
                f'string sourceGridCaseDir = "{producer_run.resolve()}";',
                config,
            )
            self.assertIn('string gridFileName = "results/grid.ofl";', config)
            self.assertIn("int maxSteps = 50;", config)
            self.assertEqual(
                original_config.read_text(encoding="utf-8"),
                'string gridFileName = "grid/old.ofl";\nint maxSteps = 50;\n',
            )

    def test_rejects_parameter_file_outside_staged_case(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "mesh").mkdir()
            (root / "solver").mkdir()
            plan = build_case_graph_plan(
                [
                    {
                        "id": "mesh",
                        "directory": "mesh",
                        "outputs": [{"id": "grid", "path": "grid.ofl"}],
                    },
                    {
                        "id": "solver",
                        "directory": "solver",
                        "depends_on": ["mesh"],
                        "inputs": [
                            {
                                "artifact": "grid",
                                "setting": "gridFileName",
                                "parameter_file": "../outside.txt",
                            }
                        ],
                    },
                ],
                root,
            )
            step = next(step for step in plan.steps if step.case_id == "solver")
            with self.assertRaisesRegex(ValueError, "stay within its producer case|stay within the staged case"):
                prepare_case_directory(
                    step,
                    root / "run",
                    {"mesh": root / "run" / "cases" / "mesh"},
                )


if __name__ == "__main__":
    unittest.main()
