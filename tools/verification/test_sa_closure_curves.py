"""Theory curve helpers for SA intermediate verification."""

from __future__ import annotations

import unittest

from sa_closure_curves import sa_fw, sa_fv1, sa_g


class SaClosureCurvesTest(unittest.TestCase):
    def test_fw_at_r_zero_is_zero(self) -> None:
        self.assertEqual(sa_fw(0.0), 0.0)

    def test_fw_at_r_one_is_one(self) -> None:
        # g(1)=1, fw(1)=1 for standard constants
        self.assertAlmostEqual(sa_fw(1.0), 1.0, places=12)

    def test_fv1_large_chi_near_one(self) -> None:
        self.assertGreater(sa_fv1(100.0), 0.999)


if __name__ == "__main__":
    unittest.main()
