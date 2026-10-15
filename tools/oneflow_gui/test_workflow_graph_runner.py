"""Tests for sequential graph execution and failure propagation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.workflow_plan import build_case_graph_plan
from tools.oneflow_gui.workflow_graph_runner import run_case_graph


class WorkflowGraphRunnerTest(unittest.TestCase):
    def _plan(self, root: Path):
        return build_case_graph_plan(
            [
                {
                    "id": "solver",
                    "directory": "solver",
                    "depends_on": ["mesh"],
                    "inputs": [{"artifact": "grid", "setting": "gridFileName", "parameter_file": "script/cfd.txt"}],
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

    def test_solver_runs_after_mesh_and_receives_runtime_overrides(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self._plan(root)
            invoked = []

            def run_step(step):
                invoked.append(step.case_id)
                if step.case_id == "mesh":
                    output = step.directory / "results/mesh.ofl"
                    output.parent.mkdir(parents=True, exist_ok=True)
                    output.write_bytes(b"mesh")
                if step.case_id == "solver":
                    self.assertEqual(
                        step.inputs[0].runtime_overrides["gridFileName"],
                        "results/mesh.ofl",
                    )
                return "succeeded", "ok"

            result = run_case_graph(plan, run_step)
            self.assertEqual(result.status, "succeeded")
            self.assertEqual(invoked, ["independent", "mesh", "solver"])

    def test_failed_mesh_blocks_solver_but_independent_case_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            invoked = []

            def run_step(step):
                invoked.append(step.case_id)
                if step.case_id == "mesh":
                    return "failed", "mesh conversion failed"
                return "succeeded", "ok"

            result = run_case_graph(self._plan(Path(temporary)), run_step)
            statuses = {item.case_id: item.status for item in result.cases}
            self.assertEqual(statuses["mesh"], "failed")
            self.assertEqual(statuses["solver"], "skipped")
            self.assertEqual(statuses["independent"], "succeeded")
            self.assertNotIn("solver", invoked)
            self.assertEqual(result.status, "failed")

    def test_missing_grid_file_blocks_solver(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            plan = self._plan(Path(temporary))
            invoked = []
            result = run_case_graph(
                plan,
                lambda step: (invoked.append(step.case_id) or "succeeded", "ok"),
            )
            statuses = {item.case_id: item.status for item in result.cases}
            self.assertEqual(statuses["mesh"], "succeeded")
            self.assertEqual(statuses["solver"], "failed")
            self.assertNotIn("solver", invoked[2:])

    def test_cancelled_case_stops_later_independent_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            invoked = []

            def run_step(step):
                invoked.append(step.case_id)
                if step.case_id == "independent":
                    return "cancelled", "user stopped"
                return "succeeded", "ok"

            result = run_case_graph(self._plan(Path(temporary)), run_step)
            self.assertEqual(invoked, ["independent"])
            self.assertTrue(all(item.status == "cancelled" for item in result.cases))
            self.assertEqual(result.status, "cancelled")

    def test_cancelled_run_does_not_launch_pending_cases(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            invoked = []
            result = run_case_graph(
                self._plan(Path(temporary)),
                lambda step: (invoked.append(step.case_id) or "succeeded", "ok"),
                is_cancelled=lambda: True,
            )
            self.assertEqual(invoked, [])
            self.assertEqual(result.status, "cancelled")
            self.assertTrue(all(item.status == "cancelled" for item in result.cases))


if __name__ == "__main__":
    unittest.main()
