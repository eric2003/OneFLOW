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

    def test_edited_expected_contract_reports_mismatch(self) -> None:
        from copy import deepcopy

        from geometry_workbench.scenarios import list_scenario_paths, load_scenario
        from geometry_workbench.pipeline_stub import diff_tables, run_pipeline

        scenario = load_scenario(list_scenario_paths(REPO)[0])
        actual = run_pipeline(scenario, backend="stub")
        edited_expected = deepcopy(scenario["expected"])
        edited_expected["l2g"][0][0] = -999

        mismatches = diff_tables(actual, edited_expected)

        self.assertTrue(mismatches)
        self.assertIn("l2g: actual=", mismatches[0])


if __name__ == "__main__":
    unittest.main()
