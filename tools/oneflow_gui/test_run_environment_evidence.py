"""Tests for immutable per-run runtime environment evidence."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from tools.oneflow_gui.run_environment_evidence import RunEnvironmentEvidence


class RunEnvironmentEvidenceTest(unittest.TestCase):
    def test_captures_json_and_finalizes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "OneFLOW.exe"
            executable.touch()
            document = {
                "schema_version": 1,
                "oneflow": {"version": "1.0.2.0"},
                "dependencies": [],
            }
            completed = subprocess.CompletedProcess(
                args=[], returncode=0, stdout=json.dumps(document), stderr=""
            )
            with patch(
                "tools.oneflow_gui.run_environment_evidence.subprocess.run",
                return_value=completed,
            ) as run:
                evidence = RunEnvironmentEvidence(
                    run_id="run-001",
                    executable=executable,
                    run_directory=root / "executions" / "run-001",
                    command_template=[str(executable), "r", "<case-directory>"],
                    working_directory=root,
                    project_directory=root / "case",
                    execution_mode="default_workspace",
                )
                original_snapshot = evidence.environment_path.read_bytes()
                evidence.finish(status="succeeded", return_code=0)

            run.assert_called_once()
            self.assertEqual(
                run.call_args.args[0],
                [str(executable.resolve()), "--environment", "--format", "json"],
            )
            self.assertEqual(
                json.loads(evidence.environment_path.read_text(encoding="utf-8"))["oneflow"]["version"],
                "1.0.2.0",
            )
            self.assertEqual(evidence.environment_path.read_bytes(), original_snapshot)
            manifest = json.loads(evidence.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "succeeded")
            self.assertEqual(manifest["return_code"], 0)
            self.assertEqual(manifest["environment_report"], "environment.json")
            self.assertEqual(manifest["environment_query_status"], "captured")
            self.assertEqual(manifest["command_template"][-1], "<case-directory>")

    def test_query_failure_is_recorded_without_blocking_run_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "OneFLOW"
            executable.touch()
            with patch(
                "tools.oneflow_gui.run_environment_evidence.subprocess.run",
                side_effect=subprocess.TimeoutExpired("OneFLOW", 15),
            ):
                evidence = RunEnvironmentEvidence(
                    run_id="run-timeout",
                    executable=executable,
                    run_directory=root / "run-timeout",
                    command_template=[str(executable), "r", "<case-directory>"],
                    working_directory=root,
                    project_directory=root / "case",
                    execution_mode="source",
                )
                evidence.finish(status="failed", return_code=1)

            report = json.loads(evidence.environment_path.read_text(encoding="utf-8"))
            manifest = json.loads(evidence.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(report["capture_status"], "failed")
            self.assertIn("timed out", report["error"])
            self.assertEqual(manifest["environment_query_status"], "failed")
            self.assertEqual(manifest["status"], "failed")


if __name__ == "__main__":
    unittest.main()
