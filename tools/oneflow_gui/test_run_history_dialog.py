"""Tests for runtime evidence links in the run history dialog."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from tools.oneflow_gui.run_history_dialog import RunHistoryDialog
from tools.oneflow_gui.workflow_repository import WorkflowRunRecord


class RunHistoryEvidenceDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])

    def test_selected_run_displays_environment_and_manifest_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            environment_path = root / "environment.json"
            manifest_path = root / "run-manifest.json"
            environment = {"schema_version": 1, "dependencies": [{"name": "CGNS"}]}
            manifest = {"schema_version": 1, "run_id": "run-001", "status": "succeeded"}
            environment_path.write_text(json.dumps(environment), encoding="utf-8")
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            record = WorkflowRunRecord(
                run_id="run-001",
                workflow_id="case-workbench",
                status="succeeded",
                started_at="2026-11-08T10:00:00+00:00",
                steps=({
                    "type": "case_queue",
                    "environment_report": str(environment_path),
                    "run_manifest": str(manifest_path),
                },),
            )
            with patch(
                "tools.oneflow_gui.run_history_dialog.WorkflowRepository"
            ) as repository_type:
                repository_type.return_value.list_runs.return_value = [record]
                dialog = RunHistoryDialog(root)
            self.addCleanup(dialog.close)

            details = dialog.details.toPlainText()
            self.assertIn("Runtime Environment Report:", details)
            self.assertIn(str(environment_path), details)
            self.assertIn('"CGNS"', details)
            self.assertIn("Run Manifest:", details)
            self.assertIn('"succeeded"', details)

    def test_missing_evidence_is_reported_without_hiding_run_details(self) -> None:
        record = WorkflowRunRecord(
            run_id="legacy-run",
            workflow_id="legacy-workflow",
            status="succeeded",
            started_at="2026-11-07T10:00:00+00:00",
        )
        with patch(
            "tools.oneflow_gui.run_history_dialog.WorkflowRepository"
        ) as repository_type:
            repository_type.return_value.list_runs.return_value = [record]
            dialog = RunHistoryDialog(Path("."))
        self.addCleanup(dialog.close)

        details = dialog.details.toPlainText()
        self.assertIn("Not linked by this run record.", details)
        self.assertIn("legacy-run", details)


if __name__ == "__main__":
    unittest.main()
