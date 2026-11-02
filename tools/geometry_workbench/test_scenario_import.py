"""Tests for GoogleTest-to-scenario bootstrap metadata."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from geometry_workbench.scenario_import import (
    scenario_id_for_test,
    scenario_skeleton,
    write_scenario_skeleton,
)


class ScenarioImportTest(unittest.TestCase):
    def setUp(self) -> None:
        self.record = {
            "framework": "GoogleTest",
            "suite": "FaceSearchBasicTest",
            "name": "RegisteredChildResolutionReturnsOnlyKnownChildren",
            "source": "tests/geometry/FaceSearchBasicTest.cpp",
            "line": 140,
            "discovery": "source",
        }

    def test_skeleton_preserves_test_identity_without_inventing_contract(self) -> None:
        scenario = scenario_skeleton(self.record)

        self.assertEqual(scenario["id"], "facesearchbasictest_registeredchildresolutionreturnsonlyknownchildren")
        self.assertEqual(scenario["test"]["source"], self.record["source"])
        self.assertEqual(scenario["provenance"]["sourceLine"], 140)
        self.assertEqual(scenario["input"], {})
        self.assertEqual(scenario["expected"], {})
        self.assertIsNone(scenario["result"])

    def test_write_skeleton_is_explicit_and_does_not_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            path = write_scenario_skeleton(self.record, Path(temp_dir))
            self.assertTrue(path.is_file())
            loaded = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(loaded["test"]["name"], self.record["name"])
            with self.assertRaises(FileExistsError):
                write_scenario_skeleton(self.record, Path(temp_dir))

    def test_non_gtest_record_is_rejected(self) -> None:
        record = dict(self.record, framework="unittest")
        with self.assertRaises(ValueError):
            scenario_id_for_test(record)


if __name__ == "__main__":
    unittest.main()
