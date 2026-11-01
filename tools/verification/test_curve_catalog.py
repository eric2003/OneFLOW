"""Catalog listing for Contract Workbench theory entries."""

from __future__ import annotations

import unittest
from pathlib import Path

from curve_catalog import list_theory_entries, resolve_repo_path, taxonomy_path


class CurveCatalogTest(unittest.TestCase):
    def test_taxonomy_and_ids(self) -> None:
        entries = list_theory_entries()
        ids = {e.id for e in entries}
        self.assertIn("sa_fw_r", ids)
        self.assertIn("sa_source_residual", ids)
        for e in entries:
            fam, gen, sec = taxonomy_path(e)
            self.assertEqual(fam, e.family)
            self.assertTrue(e.formula)
            self.assertTrue(e.references)

    def test_resolve_repo_path_portable(self) -> None:
        root = Path("/tmp/fake_repo")
        p = resolve_repo_path(root, "tools/verification/data/x.csv")
        self.assertIsNotNone(p)
        self.assertTrue(str(p).replace("\\", "/").endswith("tools/verification/data/x.csv"))


if __name__ == "__main__":
    unittest.main()
