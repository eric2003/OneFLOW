"""Analytic intermediate contract helpers."""

from __future__ import annotations

import unittest

from geometry_workbench.analytic_contracts import run_analytic, sa_fv1


class AnalyticContractsTest(unittest.TestCase):
    def test_sa_fv1_matches_closed_form(self) -> None:
        self.assertAlmostEqual(sa_fv1(0.0), 0.0)
        chi = 10.0
        cv1 = 7.1
        expected = chi**3 / (chi**3 + cv1**3)
        self.assertAlmostEqual(sa_fv1(chi, cv1), expected)

    def test_run_analytic_sa_fv1_scenario(self) -> None:
        scenario = {
            "runner": "analytic_sa_fv1",
            "input": {"cv1": 7.1, "chi": [1.0, 2.0]},
        }
        out = run_analytic(scenario)
        self.assertIsNotNone(out)
        self.assertEqual(len(out["fv1"]), 2)


if __name__ == "__main__":
    unittest.main()
