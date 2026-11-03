from __future__ import annotations

import unittest
from pathlib import Path

if __package__:
    from .curve_catalog import (
        list_theory_entries,
        resolve_repo_path,
        validate_catalog_api_references,
        validate_catalog_paths,
    )
else:
    from curve_catalog import (
        list_theory_entries,
        resolve_repo_path,
        validate_catalog_api_references,
        validate_catalog_paths,
    )


class CatalogTest(unittest.TestCase):
    def test_has_sa_sst_ns(self) -> None:
        genera = {e.genus for e in list_theory_entries()}
        self.assertIn("SA", genera)
        self.assertIn("SST", genera)
        self.assertIn("NS-LUSGS", genera)
        self.assertIn("Time integration", genera)

    def test_time_integral_entry_documents_stage_schedule_without_curve_claim(self) -> None:
        entry = next(
            item for item in list_theory_entries()
            if item.id == "time_integral_rk_stage_weight"
        )

        self.assertIsNone(entry.csv_relpath)
        self.assertIsNone(entry.plot_module)
        self.assertIn("1 / (s + 1)", entry.formula_plain)
        self.assertIn("Butcher tableau", entry.does_not_prove)
        self.assertIn("codes/solver/src/TimeIntegral.cpp", entry.source_relpaths)
        self.assertTrue(entry.examples)
        self.assertIn("rk_stage = 4", entry.examples[0])
        self.assertIn("not a complete Butcher", entry.examples[0])
    def test_representative_entries_have_worked_examples(self) -> None:
        entries = {item.id: item for item in list_theory_entries()}
        for entry_id in ("sa_fv1_chi", "sst_f1", "time_integral_rk_stage_weight"):
            self.assertTrue(entries[entry_id].examples, entry_id)

    def test_latex_and_sketch_present(self) -> None:
        for e in list_theory_entries():
            self.assertTrue(e.formula_latex)
            self.assertIn("GATED", e.code_sketch)

    def test_portable_path(self) -> None:
        p = resolve_repo_path(Path("/repo"), "tools/verification/data/x.csv")
        self.assertTrue(str(p).replace("\\", "/").endswith("tools/verification/data/x.csv"))

    def test_all_catalog_api_references_exist(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        issues = validate_catalog_api_references(repository_root)
        self.assertEqual(issues, [], "\n".join(issues))

    def test_all_catalog_references_exist(self) -> None:
        repository_root = Path(__file__).resolve().parents[2]
        issues = validate_catalog_paths(repository_root)
        self.assertEqual(issues, [], "\n".join(issues))


if __name__ == "__main__":
    unittest.main()
