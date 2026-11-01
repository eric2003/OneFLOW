from __future__ import annotations

import unittest
from pathlib import Path

from curve_catalog import (
    list_theory_entries,
    resolve_repo_path,
    validate_catalog_paths,
)


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

    def test_all_catalog_references_exist(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        issues = validate_catalog_paths(repository_root)
        self.assertEqual(issues, [], "\\n".join(issues))


if __name__ == "__main__":
    unittest.main()
