#!/usr/bin/env python3
import unittest
from pathlib import Path

from pipeline_stub import diff_tables, run_pipeline
from scenarios import load_scenario

ROOT = Path(__file__).resolve().parent


class WorkbenchStubTest(unittest.TestCase):
    def test_multi_zone_fixture_matches_declared_tables(self) -> None:
        # This legacy scenario is a display fixture, not a native algorithm run.
        scenario = load_scenario(ROOT / "scenarios" / "two_zone_parent_split.json")
        actual = run_pipeline(scenario, backend="stub")
        self.assertEqual(diff_tables(actual, scenario["expected"]), [])

    def test_registered_child_resolution_computes_from_input(self) -> None:
        scenario = load_scenario(ROOT / "scenarios" / "registered_child_resolution.json")
        actual = run_pipeline(scenario, backend="stub")
        self.assertEqual(diff_tables(actual, scenario["expected"]), [])
        self.assertEqual(
            actual["registeredChildren"],
            [
                {"globalFaceId": 0, "candidateIndex": 0},
                {"globalFaceId": 1, "candidateIndex": 1},
            ],
        )


if __name__ == "__main__":
    unittest.main()
