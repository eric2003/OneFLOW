"""Catalog listing for Contract Workbench verification curves."""

from __future__ import annotations

import unittest
from pathlib import Path

from curve_catalog import list_theory_curves, csv_path


class CurveCatalogTest(unittest.TestCase):
    def test_lists_sa_curves(self) -> None:
        ids = {c.id for c in list_theory_curves()}
        self.assertIn("sa_fw_r", ids)
        self.assertIn("sa_fv2_chi", ids)

    def test_csv_path_under_data(self) -> None:
        root = Path(__file__).resolve().parents[2]
        c = list_theory_curves()[0]
        path = csv_path(root, c)
        self.assertEqual(path.name, c.csv_name)


if __name__ == "__main__":
    unittest.main()
