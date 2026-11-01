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

    def test_verification_python_sources_are_ascii(self) -> None:
        verification_root = Path(__file__).resolve().parent
        source_paths = list(verification_root.glob("*.py"))
        source_paths.append(
            verification_root.parent
            / "oneflow_gui"
            / "verification_curve_dialog.py"
        )
        for source_path in source_paths:
            with self.subTest(path=source_path):
                try:
                    source_path.read_bytes().decode("ascii")
                except UnicodeDecodeError as exc:
                    self.fail(f"{source_path} contains non-ASCII source text: {exc}")


if __name__ == "__main__":
    unittest.main()
