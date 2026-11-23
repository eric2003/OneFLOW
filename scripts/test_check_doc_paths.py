"""Regression tests for the documentation path checker."""
from __future__ import annotations
import contextlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path

CHECKER_PATH = Path(__file__).with_name("check-doc-paths.py")
SPEC = importlib.util.spec_from_file_location("check_doc_paths", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)


class DocumentationPathCheckerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.old_root = CHECKER.REPO_ROOT
        self.old_allowlist = CHECKER.ALLOWLIST_PATH
        CHECKER.REPO_ROOT = self.root
        CHECKER.ALLOWLIST_PATH = self.root / "scripts" / "doc-path-allowlist.txt"

    def tearDown(self) -> None:
        CHECKER.REPO_ROOT = self.old_root
        CHECKER.ALLOWLIST_PATH = self.old_allowlist
        self.temp.cleanup()

    def write(self, name: str, content: str) -> None:
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def run_checker(self) -> tuple[int, str]:
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = CHECKER.main()
        return result, output.getvalue()

    def test_existing_docs_reference_passes(self) -> None:
        self.write("docs/source/mesh.rst", "Mesh architecture")
        self.write("README.md", "See docs/source/mesh.rst for details.\n")
        result, output = self.run_checker()
        self.assertEqual(result, 0, output)
        self.assertIn("OK:", output)

    def test_legacy_doc_reference_fails_once(self) -> None:
        self.write("README.md", "See doc/source/mesh.rst for details.\n")
        result, output = self.run_checker()
        self.assertEqual(result, 1)
        self.assertIn("legacy doc/ references (1)", output)
        self.assertEqual(output.count("README.md:1 doc/source/mesh.rst"), 1)

    def test_missing_docs_target_fails(self) -> None:
        self.write("README.md", "See docs/source/missing.rst.\n")
        result, output = self.run_checker()
        self.assertEqual(result, 1)
        self.assertIn("missing docs/ references (1)", output)

    def test_checker_source_is_not_scanned_against_its_own_patterns(self) -> None:
        self.write("scripts/check-doc-paths.py", "Example: doc/source/old.rst\n")
        self.write("README.md", "No path references here.\n")
        result, output = self.run_checker()
        self.assertEqual(result, 0, output)


if __name__ == "__main__":
    unittest.main()
