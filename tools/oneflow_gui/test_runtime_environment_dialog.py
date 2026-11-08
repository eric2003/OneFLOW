"""Tests for the OneFLOW runtime environment query and dialog service."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.runtime_environment_dialog import query_runtime_environment


class RuntimeEnvironmentQueryTest(unittest.TestCase):
    def _executable(self, directory: str) -> Path:
        executable = Path(directory) / "OneFLOW"
        executable.write_text("test executable placeholder", encoding="ascii")
        return executable

    def _document(self) -> dict:
        return {
            "schema_version": 1,
            "oneflow": {"version": "1.2.3"},
            "build": {
                "compiler": "MSVC",
                "compiler_version": "19.0",
                "configuration": "Debug",
            },
            "dependencies": [
                {
                    "name": "CGNS",
                    "enabled": True,
                    "status": "available",
                    "runtime_version": {
                        "value": None,
                        "source": "unavailable",
                    },
                    "header_version": {
                        "value": "4.5.2",
                        "source": "CGNS_VERSION_compile_time_macro",
                    },
                    "version_compatibility": {
                        "status": "assumed_match",
                        "basis": "header_version",
                        "detail": "Runtime version was not independently obtained.",
                    },
                    "runtime_library": {
                        "module_path": "C:/dev/cgns/bin/cgnsdll.dll",
                        "linkage": "shared",
                    },
                },
                {
                    "name": "HDF5",
                    "enabled": False,
                    "status": "not_detected",
                    "runtime_version": {
                        "value": None,
                        "source": "unavailable",
                    },
                    "runtime_version_status": "runtime_api_symbol_not_found",
                    "header_version": None,
                    "version_compatibility": {
                        "status": "unknown",
                        "basis": "insufficient_version_evidence",
                        "detail": "A runtime version alone does not establish the header version.",
                    },
                    "runtime_library": {
                        "module_path": "",
                        "linkage": "unavailable",
                    },
                    "detail": "H5get_libversion was not exported by a loaded process module.",
                }
            ],
        }

    def test_invokes_oneflow_and_returns_validated_environment_document(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = self._executable(directory)
            document = self._document()
            completed = subprocess.CompletedProcess(
                args=[],
                returncode=0,
                stdout=json.dumps(document),
                stderr="",
            )
            with patch(
                "tools.oneflow_gui.runtime_environment_dialog.subprocess.run",
                return_value=completed,
            ) as run:
                result = query_runtime_environment(executable)

            self.assertEqual(result, document)
            self.assertEqual(
                result["dependencies"][0]["version_compatibility"]["status"],
                "assumed_match",
            )
            self.assertEqual(
                result["dependencies"][1]["version_compatibility"]["status"],
                "unknown",
            )
            self.assertIsNone(result["dependencies"][1]["header_version"])
            arguments, keywords = run.call_args
            self.assertEqual(
                arguments[0],
                [
                    str(executable.resolve()),
                    "--environment",
                    "--format",
                    "json",
                ],
            )
            self.assertEqual(keywords["cwd"], executable.resolve().parent)
            self.assertEqual(keywords["timeout"], 15.0)
            self.assertFalse(keywords["check"])

    def test_reports_nonzero_exit_with_process_diagnostic(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = self._executable(directory)
            completed = subprocess.CompletedProcess(
                args=[],
                returncode=3,
                stdout="",
                stderr="CGNS runtime query failed",
            )
            with patch(
                "tools.oneflow_gui.runtime_environment_dialog.subprocess.run",
                return_value=completed,
            ):
                with self.assertRaisesRegex(
                    RuntimeError, "CGNS runtime query failed"
                ):
                    query_runtime_environment(executable)

    def test_rejects_invalid_json(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = self._executable(directory)
            completed = subprocess.CompletedProcess(
                args=[], returncode=0, stdout="{", stderr=""
            )
            with patch(
                "tools.oneflow_gui.runtime_environment_dialog.subprocess.run",
                return_value=completed,
            ):
                with self.assertRaisesRegex(ValueError, "invalid environment JSON"):
                    query_runtime_environment(executable)

    def test_rejects_unknown_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = self._executable(directory)
            document = self._document()
            document["schema_version"] = 2
            completed = subprocess.CompletedProcess(
                args=[],
                returncode=0,
                stdout=json.dumps(document),
                stderr="",
            )
            with patch(
                "tools.oneflow_gui.runtime_environment_dialog.subprocess.run",
                return_value=completed,
            ):
                with self.assertRaisesRegex(ValueError, "Unsupported"):
                    query_runtime_environment(executable)

    def test_reports_missing_executable_before_process_launch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "missing-OneFLOW"
            with patch(
                "tools.oneflow_gui.runtime_environment_dialog.subprocess.run"
            ) as run:
                with self.assertRaises(FileNotFoundError):
                    query_runtime_environment(missing)
            run.assert_not_called()

    def test_converts_process_timeout_to_actionable_error(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            executable = self._executable(directory)
            timeout_error = subprocess.TimeoutExpired(
                cmd=[str(executable)], timeout=2.0
            )
            with patch(
                "tools.oneflow_gui.runtime_environment_dialog.subprocess.run",
                side_effect=timeout_error,
            ):
                with self.assertRaisesRegex(RuntimeError, "exceeded 15 seconds"):
                    query_runtime_environment(executable)


if __name__ == "__main__":
    unittest.main()
