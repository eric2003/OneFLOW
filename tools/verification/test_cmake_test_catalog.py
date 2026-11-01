from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from .cmake_test_catalog import discover_cmake_test_targets


class CMakeTestCatalogTest(unittest.TestCase):
    def test_maps_source_to_target_and_ctest_discovery_declaration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tests = root / "tests" / "runtime"
            tests.mkdir(parents=True)
            source = tests / "lifecycle_test.cpp"
            source.write_text("TEST(Lifecycle, ReleasesState) {}\n", encoding="utf-8")
            (tests / "CMakeLists.txt").write_text(
                "add_executable(lifecycle_test lifecycle_test.cpp)\n"
                "gtest_discover_tests(lifecycle_test DISCOVERY_MODE PRE_TEST)\n",
                encoding="utf-8",
            )

            records = discover_cmake_test_targets(
                root, "tests/runtime/lifecycle_test.cpp"
            )

            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].target, "lifecycle_test")
            self.assertEqual(records[0].cmake_file, "tests/runtime/CMakeLists.txt")
            self.assertTrue(records[0].discovery_registered)
            self.assertEqual(records[0].discovery_mode, "PRE_TEST")

    def test_does_not_infer_discovery_from_executable_declaration(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tests = root / "tests"
            tests.mkdir()
            (tests / "sample.cpp").write_text("TEST(Sample, Runs) {}\n", encoding="utf-8")
            (tests / "CMakeLists.txt").write_text(
                "add_executable(sample_test sample.cpp)\n",
                encoding="utf-8",
            )

            records = discover_cmake_test_targets(root, "tests/sample.cpp")

            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].target, "sample_test")
            self.assertFalse(records[0].discovery_registered)
            self.assertEqual(records[0].discovery_mode, "")

    def test_does_not_match_unrelated_or_variable_indirect_source(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            tests = root / "tests"
            tests.mkdir()
            (tests / "sample.cpp").write_text("", encoding="utf-8")
            (tests / "other.cpp").write_text("", encoding="utf-8")
            (tests / "CMakeLists.txt").write_text(
                "set(TEST_SOURCES sample.cpp)\n"
                "add_executable(sample_test ${TEST_SOURCES})\n",
                encoding="utf-8",
            )

            self.assertEqual(
                discover_cmake_test_targets(root, "tests/sample.cpp"),
                (),
            )
            self.assertEqual(
                discover_cmake_test_targets(root, "tests/other.cpp"),
                (),
            )

    def test_rejects_source_outside_repository(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.assertEqual(
                discover_cmake_test_targets(root, "../outside.cpp"),
                (),
            )


if __name__ == "__main__":
    unittest.main()
