from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.verification.preflight import (
    _changed_files,
    _cpp_test_target_diagnostic,
    _stale_cpp_evidence_diagnostic,
    find_cpp_test_candidates,
    find_python_test_candidates,
    format_report,
)


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
            subprocess.run(("git", "config", "user.email", "preflight@example.invalid"), cwd=root, check=True)
            subprocess.run(("git", "config", "user.name", "Preflight Test"), cwd=root, check=True)
            (root / "tracked.txt").write_text("base", encoding="utf-8")
            subprocess.run(("git", "add", "tracked.txt"), cwd=root, check=True)
            subprocess.run(("git", "commit", "-qm", "initial"), cwd=root, check=True)
            (root / "new.txt").write_text("new", encoding="utf-8")
            self.assertEqual(_changed_files(root), ("new.txt",))



    def test_stale_python_test_path_finds_exact_declaration_after_move(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            moved = root / "grid" / "tests" / "python" / "test_moved.py"
            moved.parent.mkdir(parents=True)
            moved.write_text(
                "class Stl2DGmshTest(unittest.TestCase):\n"
                "    def test_outer_domain_allows_shared_physical_group_for_multiple_sides(self):\n"
                "        pass\n",
                encoding="utf-8",
            )
            (root / "build").mkdir()
            (root / "build" / "ignored.py").write_text(
                "class Stl2DGmshTest:\n"
                "    def test_outer_domain_allows_shared_physical_group_for_multiple_sides(self): pass\n",
                encoding="utf-8",
            )
            self.assertEqual(
                find_python_test_candidates(
                    root,
                    "Stl2DGmshTest.test_outer_domain_allows_shared_physical_group_for_multiple_sides",
                ),
                ("grid/tests/python/test_moved.py",),
            )


    def test_stale_google_test_path_finds_declaration_after_move(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            moved = root / "tests" / "geometry" / "StlMeshIOTest.cpp"
            moved.parent.mkdir(parents=True)
            moved.write_text(
                "// TEST(OldSuite, MovedCase)\n"
                "TEST_F(StlMeshIOTest, PreservesBoundaryIdentity) { }\n"
                "const char* example = \"TEST(StlMeshIOTest, PreservesBoundaryIdentity)\";\n",
                encoding="utf-8",
            )
            build = root / "build"
            build.mkdir()
            (build / "generated.cpp").write_text(
                "TEST_F(StlMeshIOTest, PreservesBoundaryIdentity) { }\n",
                encoding="utf-8",
            )
            self.assertEqual(
                find_cpp_test_candidates(
                    root, "StlMeshIOTest.PreservesBoundaryIdentity"
                ),
                ("tests/geometry/StlMeshIOTest.cpp",),
            )
            diagnostic = _stale_cpp_evidence_diagnostic(
                root,
                "StlMeshIOTest.PreservesBoundaryIdentity",
                "missing repository file: tests/old/StlMeshIOTest.cpp",
            )
            self.assertIn(
                "matching GoogleTest declaration found at "
                "tests/geometry/StlMeshIOTest.cpp",
                diagnostic,
            )

    def test_google_test_target_without_ctest_discovery_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tests = root / "tests" / "geometry"
            tests.mkdir(parents=True)
            (tests / "StlMeshIOTest.cpp").write_text("TEST(StlMeshIOTest, ReadsSurface) {}", encoding="utf-8")
            (tests / "CMakeLists.txt").write_text("add_executable(StlMeshIOTest StlMeshIOTest.cpp)", encoding="utf-8")
            diagnostic = _cpp_test_target_diagnostic(root, "tests/geometry/StlMeshIOTest.cpp")
            self.assertIn("StlMeshIOTest", diagnostic)
            self.assertIn("gtest_discover_tests", diagnostic)

    def test_google_test_target_with_discovery_is_consistent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tests = root / "tests" / "geometry"
            tests.mkdir(parents=True)
            (tests / "StlMeshIOTest.cpp").write_text("TEST(StlMeshIOTest, ReadsSurface) {}", encoding="utf-8")
            (tests / "CMakeLists.txt").write_text(
                "add_executable(StlMeshIOTest StlMeshIOTest.cpp) "
                "gtest_discover_tests(StlMeshIOTest DISCOVERY_MODE PRE_TEST)",
                encoding="utf-8",
            )
            self.assertEqual(_cpp_test_target_diagnostic(root, "tests/geometry/StlMeshIOTest.cpp"), "")

    def test_stale_google_test_path_reports_ambiguous_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for folder in ("one", "two"):
                path = root / folder / "test_case.cc"
                path.parent.mkdir(parents=True)
                path.write_text(
                    "TEST(MovedSuite, MovedCase) { }\n", encoding="utf-8"
                )
            self.assertEqual(
                find_cpp_test_candidates(root, "MovedSuite.MovedCase"),
                ("one/test_case.cc", "two/test_case.cc"),
            )

    def test_stale_python_test_path_reports_all_ambiguous_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            identity = "MovedTest.test_identity"
            for folder in ("one", "two"):
                path = root / folder / "test_case.py"
                path.parent.mkdir(parents=True)
                path.write_text(
                    "class MovedTest:\n"
                    "    def test_identity(self): pass\n",
                    encoding="utf-8",
                )
            self.assertEqual(
                find_python_test_candidates(root, identity),
                ("one/test_case.py", "two/test_case.py"),
            )


if __name__ == "__main__":
    unittest.main()
