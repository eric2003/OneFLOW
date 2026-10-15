"""Tests for isolated sequential graph execution and failure propagation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_graph_runner import run_case_graph
from tools.oneflow_gui.workflow_plan import build_case_graph_plan


class WorkflowGraphRunnerTest(unittest.TestCase):
    def _plan(self, root: Path):
        (root / "solver" / "script").mkdir(parents=True)
        (root / "independent").mkdir()
        (root / "mesh").mkdir()
        (root / "solver" / "script" / "cfd.txt").write_text(
            'string gridFileName = "grid/old.ofl";\n', encoding="utf-8"
        )
        return build_case_graph_plan(
            [
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
                {"id": "independent", "directory": "independent"},
                {
                    "id": "mesh",
                    "directory": "mesh",
                    "outputs": [{"id": "grid", "path": "results/mesh.ofl"}],
                },
            ],
            root,
        )

    def test_solver_runs_after_mesh_with_staged_artifact_binding(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self._plan(root)
            run_root = root / "run"
            invoked = []

            def run_step(step):
                invoked.append(step.case_id)
                self.assertEqual(step.directory, run_root / "cases" / step.case_id)
                if step.case_id == "mesh":
                    output = step.directory / "results" / "mesh.ofl"
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_bytes(b"mesh")
                if step.case_id == "solver":
                    binding = step.inputs[0]
                    producer_directory = (run_root / "cases" / "mesh").resolve()
                    self.assertEqual(binding.source_case_directory, producer_directory)
                    self.assertTrue(binding.source_path.is_file())
                    self.assertEqual(
                        binding.runtime_overrides["sourceGridCaseDir"],
                        str(producer_directory),
                    )
                    self.assertEqual(
                        binding.runtime_overrides["gridFileName"],
                        "results/mesh.ofl",
                    )
                    config = (step.directory / "script" / "cfd.txt").read_text(encoding="utf-8")
                    escaped_producer_path = str(producer_directory).replace(chr(92), chr(92) * 2)
                    self.assertIn(
                        f'string sourceGridCaseDir = "{escaped_producer_path}";',
                        config,
                    )
                    self.assertIn('string gridFileName = "results/mesh.ofl";', config)
                return "succeeded", "ok"

            result = run_case_graph(plan, run_step, run_root=run_root)
            self.assertEqual(result.status, "succeeded")
            self.assertEqual(invoked, ["independent", "mesh", "solver"])
            self.assertEqual(
                (root / "solver" / "script" / "cfd.txt").read_text(encoding="utf-8"),
                'string gridFileName = "grid/old.ofl";\n',
            )

    def test_failed_mesh_blocks_solver_but_independent_case_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self._plan(root)
            invoked = []

            def run_step(step):
                invoked.append(step.case_id)
                if step.case_id == "mesh":
                    return "failed", "mesh conversion failed"
                return "succeeded", "ok"

            result = run_case_graph(plan, run_step, run_root=root / "run")
            statuses = {item.case_id: item.status for item in result.cases}
            self.assertEqual(statuses["mesh"], "failed")
            self.assertEqual(statuses["solver"], "skipped")
            self.assertEqual(statuses["independent"], "succeeded")
            self.assertNotIn("solver", invoked)
            self.assertEqual(result.status, "failed")

    def test_missing_grid_file_blocks_solver(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self._plan(root)
            invoked = []
            result = run_case_graph(
                plan,
                lambda step: (invoked.append(step.case_id) or "succeeded", "ok"),
                run_root=root / "run",
            )
            statuses = {item.case_id: item.status for item in result.cases}
            self.assertEqual(statuses["mesh"], "succeeded")
            self.assertEqual(statuses["solver"], "failed")
            self.assertNotIn("solver", invoked)

    def test_cancelled_case_stops_later_independent_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self._plan(root)
            invoked = []

            def run_step(step):
                invoked.append(step.case_id)
                if step.case_id == "independent":
                    return "cancelled", "user stopped"
                return "succeeded", "ok"

            result = run_case_graph(plan, run_step, run_root=root / "run")
            self.assertEqual(invoked, ["independent"])
            self.assertTrue(all(item.status == "cancelled" for item in result.cases))
            self.assertEqual(result.status, "cancelled")

    def test_cancelled_run_does_not_stage_or_launch_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self._plan(root)
            invoked = []
            run_root = root / "run"
            result = run_case_graph(
                plan,
                lambda step: (invoked.append(step.case_id) or "succeeded", "ok"),
                run_root=run_root,
                is_cancelled=lambda: True,
            )
            self.assertEqual(invoked, [])
            self.assertEqual(result.status, "cancelled")
            self.assertTrue(all(item.status == "cancelled" for item in result.cases))
            self.assertFalse(run_root.exists())


if __name__ == "__main__":
    unittest.main()
