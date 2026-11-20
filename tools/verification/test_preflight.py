from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.verification.preflight import format_report, _changed_files


class ChangePreflightTest(unittest.TestCase):
    def test_report_lists_each_consistency_failure(self):
        report = format_report(
            ("missing contract source", "unregistered Python test"),
            ("tools/verification/preflight.py",),
        )
        self.assertIn("Changed files: 1", report)
        self.assertIn("FAIL: 2 issue(s)", report)
        self.assertIn("missing contract source", report)
        self.assertIn("unregistered Python test", report)

    def test_report_has_explicit_pass_state(self):
        self.assertIn("PASS:", format_report(()))

    def test_changed_file_discovery_includes_untracked_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            import subprocess
            subprocess.run(("git", "init", "-q"), cwd=root, check=True)
            (root / "new.txt").write_text("new", encoding="utf-8")
            self.assertEqual(_changed_files(root), ("new.txt",))


if __name__ == "__main__":
    unittest.main()
