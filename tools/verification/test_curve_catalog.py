from __future__ import annotations

import unittest
from pathlib import Path

from curve_catalog import list_theory_entries, resolve_repo_path


class CatalogTest(unittest.TestCase):
    def test_has_sa_sst_ns(self) -> None:
        genera = {e.genus for e in list_theory_entries()}
        self.assertIn("SA", genera)
        self.assertIn("SST", genera)
        self.assertIn("NS-LUSGS", genera)

    def test_latex_and_sketch_present(self) -> None:
        for e in list_theory_entries():
            self.assertTrue(e.formula_latex)
            self.assertIn("GATED", e.code_sketch)

    def test_portable_path(self) -> None:
        p = resolve_repo_path(Path("/repo"), "tools/verification/data/x.csv")
        self.assertTrue(str(p).replace("\\", "/").endswith("tools/verification/data/x.csv"))


if __name__ == "__main__":
    unittest.main()
