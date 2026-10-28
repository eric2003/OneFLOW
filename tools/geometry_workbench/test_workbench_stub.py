#!/usr/bin/env python3
import unittest
from pathlib import Path

from pipeline_stub import diff_tables, run_pipeline
from scenarios import load_scenario

ROOT = Path(__file__).resolve().parent


class WorkbenchStubTest(unittest.TestCase):
    def test_golden_scenario_matches_itself(self) -> None:
        scenario = load_scenario(ROOT / "scenarios" / "two_zone_parent_split.json")
        actual = run_pipeline(scenario, backend="stub")
        self.assertEqual(diff_tables(actual, scenario["expected"]), [])


if __name__ == "__main__":
    unittest.main()
