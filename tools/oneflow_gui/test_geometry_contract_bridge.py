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
        scenario = load_scenario(next(path for path in paths if path.stem == "two_zone_parent_split"))
        actual = run_pipeline(scenario, backend="stub")
        self.assertEqual(diff_tables(actual, scenario["expected"]), [])

    def test_stub_output_is_fixture_backed_not_input_driven(self) -> None:
        from copy import deepcopy

        from geometry_workbench.scenarios import list_scenario_paths, load_scenario
        from geometry_workbench.pipeline_stub import run_pipeline

        scenario = load_scenario(next(path for path in list_scenario_paths(REPO) if path.stem == "two_zone_parent_split"))
        baseline = run_pipeline(scenario, backend="stub")
        edited = deepcopy(scenario)
        edited["input"]["l2g"][0][0] = 9876

        self.assertEqual(run_pipeline(edited, backend="stub"), baseline)

    def test_edited_expected_contract_reports_mismatch(self) -> None:
        from copy import deepcopy

        from geometry_workbench.scenarios import list_scenario_paths, load_scenario
        from geometry_workbench.pipeline_stub import diff_tables, run_pipeline

        scenario = load_scenario(next(path for path in list_scenario_paths(REPO) if path.stem == "two_zone_parent_split"))
        actual = run_pipeline(scenario, backend="stub")
        edited_expected = deepcopy(scenario["expected"])
        edited_expected["l2g"][0][0] = -999

        mismatches = diff_tables(actual, edited_expected)

        self.assertTrue(mismatches)
        self.assertIn("l2g: actual=", mismatches[0])

    def test_partial_child_match_scenario_points_to_cpp_regression(self) -> None:
        from geometry_workbench.scenarios import list_scenario_paths, load_scenario

        scenarios = [
            load_scenario(path)
            for path in list_scenario_paths(REPO)
            if path.stem == "registered_child_partial_match"
        ]
        self.assertEqual(len(scenarios), 1)
        scenario = scenarios[0]
        self.assertEqual(
            scenario["test"]["name"],
            "RegisteredChildResolutionPreservesCandidateIndexForPartialMatch",
        )
        self.assertEqual(
            scenario["integrationTest"]["name"],
            "QuadSplitUsesOriginalCandidateIndexForPartialMatch",
        )
        self.assertEqual(scenario["expected"]["childFaceIds"], [0])
        self.assertEqual(scenario["expected"]["candidateIndices"], [1])
        self.assertIn("fixture backend", " ".join(scenario["notes"]))


if __name__ == "__main__":
    unittest.main()
