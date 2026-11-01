"""Tests for structured test failure diagnostics."""

from __future__ import annotations

import unittest

from tools.oneflow_gui.test_diagnostics import parse_test_diagnostics


class TestDiagnosticParserTest(unittest.TestCase):
    def test_unittest_failure_captures_case_location_exception_and_context(self) -> None:
        records = parse_test_diagnostics(
            [
                "FAIL: test_contract_is_verified (sample.TestContract)",
                "Traceback (most recent call last):",
                '  File "tests/test_contract.py", line 42, in test_contract_is_verified',
                "    self.assertTrue(result)",
                "AssertionError: expected true",
                "Ran 1 test in 0.01s",
                "FAILED (failures=1)",
            ],
            task_label="Contract catalog tests",
            framework="unittest",
            configuration="Debug",
            working_directory="D:/repo",
            started_at="2026-11-02T10:00:00+08:00",
            finished_at="2026-11-02T10:00:01+08:00",
            exit_code=1,
            contract_ids=("runtime/process-lifecycle",),
        )

        self.assertEqual(len(records), 1)
        item = records[0]
        self.assertEqual(item.test_id, "test_contract_is_verified (sample.TestContract)")
        self.assertEqual(item.source_file, "tests/test_contract.py")
        self.assertEqual(item.source_line, 42)
        self.assertEqual(item.exception_type, "AssertionError")
        self.assertEqual(item.message, "expected true")
        self.assertEqual(item.configuration, "Debug")
        self.assertEqual(item.working_directory, "D:/repo")
        self.assertEqual(item.exit_code, 1)
        self.assertEqual(item.contract_ids, ("runtime/process-lifecycle",))
        self.assertEqual(item.verification_impact, "requires_review")
        self.assertEqual(item.to_dict()["contract_ids"], ["runtime/process-lifecycle"])

    def test_unittest_error_without_source_location_keeps_unknown_fields_empty(self) -> None:
        records = parse_test_diagnostics(
            ["ERROR: test_import (sample.TestImport)", "ImportError: missing module"],
            task_label="Python tests",
            framework="unittest",
        )
        self.assertEqual(len(records), 1)
        self.assertIsNone(records[0].source_file)
        self.assertIsNone(records[0].source_line)
        self.assertEqual(records[0].exception_type, "ImportError")
        self.assertEqual(records[0].message, "missing module")
        self.assertEqual(records[0].verification_impact, "unlinked")

    def test_googletest_failure_captures_cpp_location(self) -> None:
        records = parse_test_diagnostics(
            [
                "[ RUN      ] FaceTopologyTest.RejectsInvalidFace",
                "D:/repo/tests/geometry/FaceTopologyTest.cpp:87: Failure",
                "Expected: face.IsValid()",
                "[  FAILED  ] FaceTopologyTest.RejectsInvalidFace (2 ms)",
            ],
            task_label="C++ unit tests (CTest)",
            framework="GoogleTest",
            configuration="Release",
            working_directory="D:/repo/build",
            exit_code=8,
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].test_name, "FaceTopologyTest.RejectsInvalidFace")
        self.assertEqual(records[0].source_file, "D:/repo/tests/geometry/FaceTopologyTest.cpp")
        self.assertEqual(records[0].source_line, 87)
        self.assertEqual(records[0].configuration, "Release")

    def test_passed_output_does_not_create_failure_diagnostic(self) -> None:
        records = parse_test_diagnostics(
            ["[ RUN      ] SampleTest.Passes", "[       OK ] SampleTest.Passes (1 ms)"],
            task_label="C++ unit tests",
            framework="GoogleTest",
            exit_code=0,
        )
        self.assertEqual(records, ())

    def test_contract_ids_are_deduplicated_and_explicit(self) -> None:
        records = parse_test_diagnostics(
            ["FAIL: test_x (sample.Test)", "AssertionError: no"],
            task_label="task",
            framework="unittest",
            contract_ids=("case/plate", "case/plate"),
        )
        self.assertEqual(records[0].contract_ids, ("case/plate",))


if __name__ == "__main__":
    unittest.main()
