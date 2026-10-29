"""Persistence tests for editable geometry scenario contracts."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
import sys

TOOLS = Path(__file__).resolve().parents[1]
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from geometry_workbench.scenarios import save_scenario_contract


class ScenarioContractPersistenceTest(unittest.TestCase):
    def test_save_updates_tables_and_preserves_provenance_and_result(self) -> None:
        scenario = {
            "id": "sample",
            "test": {"framework": "GoogleTest", "suite": "Suite", "name": "Case"},
            "provenance": {"kind": "gtest-inventory-skeleton", "sourceLine": 12},
            "input": {"old": 1},
            "expected": {"old": 2},
            "result": None,
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "sample.json"
            path.write_text(json.dumps(scenario), encoding="utf-8")
            save_scenario_contract(path, scenario, {"nodes": [1, 2]}, {"count": 2})
            saved = json.loads(path.read_text(encoding="utf-8"))

        self.assertEqual(saved["input"], {"nodes": [1, 2]})
        self.assertEqual(saved["expected"], {"count": 2})
        self.assertEqual(saved["test"], scenario["test"])
        self.assertEqual(saved["provenance"], scenario["provenance"])
        self.assertIsNone(saved["result"])

    def test_save_rejects_non_object_tables(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "sample.json"
            with self.assertRaises(ValueError):
                save_scenario_contract(path, {}, [], {})


if __name__ == "__main__":
    unittest.main()
