"""Bridge tests: Workbench can load geometry_workbench scenarios without Qt widgets."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))


class GeometryContractBridgeTest(unittest.TestCase):
    def test_list_and_run_stub_scenario(self) -> None:
        from geometry_workbench.scenarios import list_scenario_paths, load_scenario
        from geometry_workbench.pipeline_stub import diff_tables, run_pipeline

        paths = list_scenario_paths(REPO)
        self.assertTrue(paths, "expected at least one scenario JSON")
        scenario = load_scenario(paths[0])
        actual = run_pipeline(scenario, backend="stub")
        self.assertEqual(diff_tables(actual, scenario["expected"]), [])


if __name__ == "__main__":
    unittest.main()
