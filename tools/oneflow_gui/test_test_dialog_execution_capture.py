"""Integration tests for execution evidence at the test worker boundary."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication

from tools.oneflow_gui.case_repository import load_cases, select_cases_by_execution_kind
from tools.oneflow_gui.test_dialog import (
    CaseDatabaseSelectionDialog,
    TestCommandWorker,
)
from tools.oneflow_gui.test_runner import TestCommand


class DevelopmentCaseSelectorModeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])

    def _visible_case_ids(self, execution_kind: str) -> set[str]:
        all_cases = load_cases(Path(__file__).resolve().parents[2])
        active_cases = select_cases_by_execution_kind(all_cases, execution_kind)
        dialog = CaseDatabaseSelectionDialog(active_cases)
        visible_ids = set()
        for group_index in range(dialog.case_tree.topLevelItemCount()):
            group = dialog.case_tree.topLevelItem(group_index)
            for case_index in range(group.childCount()):
                visible_ids.add(
                    group.child(case_index).data(0, Qt.ItemDataRole.UserRole)
                )
        dialog.close()
        return visible_ids

    def test_native_mode_selector_never_displays_python_cases(self) -> None:
        all_cases = load_cases(Path(__file__).resolve().parents[2])
        native_ids = {
            case.case_id for case in select_cases_by_execution_kind(all_cases, "oneflow")
        }
        python_ids = {
            case.case_id for case in select_cases_by_execution_kind(all_cases, "python")
        }
        self.assertEqual(self._visible_case_ids("oneflow"), native_ids)
        self.assertFalse(self._visible_case_ids("oneflow") & python_ids)

    def test_python_mode_selector_never_displays_native_cases(self) -> None:
        all_cases = load_cases(Path(__file__).resolve().parents[2])
        native_ids = {
            case.case_id for case in select_cases_by_execution_kind(all_cases, "oneflow")
        }
        python_ids = {
            case.case_id for case in select_cases_by_execution_kind(all_cases, "python")
        }
        self.assertEqual(self._visible_case_ids("python"), python_ids)
        self.assertFalse(self._visible_case_ids("python") & native_ids)


class TestCommandWorkerEvidenceTest(unittest.TestCase):
    def _command(self, root: Path) -> TestCommand:
        return TestCommand(
            "Sample unittest suite",
            ("python", "-m", "unittest", "tools.sample_tests"),
            root,
        )

    def _process(self, return_code: int = 0) -> MagicMock:
        process = MagicMock()
        process.stdout = iter(("Ran 1 test in 0.01s\n", "OK\n"))
        process.wait.return_value = return_code
        process.__enter__.return_value = process
        process.__exit__.return_value = False
        return process

    def test_worker_persists_evidence_after_unittest_process_completes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            process = self._process()
            paths = SimpleNamespace(
                verification_execution_records=root / ".oneflow" / "verification" / "execution-records.jsonl",
                test_run_logs=root / ".oneflow" / "test-runs" / "logs",
            )
            with (
                patch(
                    "tools.oneflow_gui.test_dialog.subprocess.Popen",
                    return_value=process,
                ),
                patch(
                    "tools.oneflow_gui.test_dialog.UserDataPaths.discover",
                    return_value=paths,
                ),
                patch(
                    "tools.oneflow_gui.test_dialog.capture_python_unittest_execution",
                    return_value=SimpleNamespace(
                        record_id="record-1", outcome="passed"
                    ),
                ) as capture,
            ):
                worker = TestCommandWorker([self._command(root)], build_config="Debug")
                worker.run()

            capture.assert_called_once()
            arguments = capture.call_args.kwargs
            self.assertEqual(arguments["exit_code"], 0)
            self.assertEqual(arguments["output_lines"], [
                "Ran 1 test in 0.01s",
                "OK",
            ])
            self.assertEqual(arguments["configuration"], "Debug")
            self.assertEqual(arguments["artifact_directory"], paths.test_run_logs)
            self.assertFalse(arguments["cancelled"])

    def test_process_launch_failure_does_not_create_execution_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch(
                    "tools.oneflow_gui.test_dialog.subprocess.Popen",
                    side_effect=OSError("python executable unavailable"),
                ),
                patch(
                    "tools.oneflow_gui.test_dialog.capture_python_unittest_execution"
                ) as capture,
            ):
                worker = TestCommandWorker([self._command(root)])
                worker.run()

            capture.assert_not_called()

    def test_cancellation_before_launch_does_not_create_execution_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with (
                patch("tools.oneflow_gui.test_dialog.subprocess.Popen") as popen,
                patch(
                    "tools.oneflow_gui.test_dialog.capture_python_unittest_execution"
                ) as capture,
            ):
                worker = TestCommandWorker([self._command(root)])
                worker.request_stop()
                worker.run()

            popen.assert_not_called()
            capture.assert_not_called()


if __name__ == "__main__":
    unittest.main()
