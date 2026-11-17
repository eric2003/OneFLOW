"""Tests for declarative Case profile loading and resolution."""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.case_profile import (
    load_case_profile,
    resolve_case_profile,
)
from tools.oneflow_gui.parameter_document import ParameterDocument


class CaseProfileTest(unittest.TestCase):
    def test_flat_plate_profiles_describe_regression_and_production_intent(self):
        root = Path(__file__).resolve().parent / "case_profiles"
        regression = load_case_profile(root / "laminar-flat-plate-regression-50.json")
        production = load_case_profile(root / "laminar-flat-plate-production-35000.json")

        self.assertEqual(regression.source_case, "test/cfd/plateuns2dslau2")
        self.assertEqual(regression.run_purpose, "regression")
        self.assertEqual(production.run_purpose, "production")
        self.assertEqual(regression.maturity, "proposed")
        self.assertEqual(production.maturity, "proposed")
        self.assertEqual([(item.key, item.value) for item in regression.overrides], [("maxSteps", 50)])
        self.assertEqual([(item.key, item.value) for item in production.overrides], [("maxSteps", 35000)])

    def test_profile_resolves_against_logical_key_not_filename(self):
        profile_path = (
            Path(__file__).resolve().parent
            / "case_profiles"
            / "laminar-flat-plate-production-35000.json"
        )
        profile = load_case_profile(profile_path)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "iteration-policy.txt"
            path.write_text(
                "int maxSteps = 50; // total steps\n",
                encoding="utf-8",
            )
            document = ParameterDocument.load(path)
            plan = resolve_case_profile(profile, [document])
            self.assertEqual(plan.changes[0].requested_value, 35000)
            self.assertIn("maxSteps: 50 -> 35000", plan.preview()[0])
            self.assertEqual(document.get("maxSteps").value, 50)
            plan.apply()
            self.assertEqual(document.get("maxSteps").value, 35000)

    def test_rejects_duplicate_or_non_scalar_overrides(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.json"
            base = {
                "schema_version": 1,
                "profile_id": "bad-profile",
                "title": "Bad",
                "source_case": "test/cfd/example",
                "run_purpose": "test",
                "maturity": "proposed",
                "description": "Invalid profile fixture",
                "limitations": [],
            }
            for overrides in (
                [{"key": "maxSteps", "value": 50}, {"key": "MAXSTEPS", "value": 100}],
                [{"key": "maxSteps", "value": {"steps": 50}}],
            ):
                path.write_text(json.dumps({**base, "overrides": overrides}), encoding="utf-8")
                with self.assertRaises(ValueError):
                    load_case_profile(path)


if __name__ == "__main__":
    unittest.main()
