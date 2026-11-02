"""Tests for the unified Contract Workbench domain navigation."""

from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication

from tools.verification.execution_records import (
    ExecutionArtifact,
    ExecutionRecord,
    ExecutionRecordStore,
)

from tools.oneflow_gui.contract_workbench_dialog import (
    ContractDomainDialog,
    ContractWorkbenchDialog,
)


class ContractWorkbenchDialogTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.application = QApplication.instance() or QApplication([])
        cls.repository_root = Path(__file__).resolve().parents[2]

    def test_domain_hub_exposes_current_contract_domains(self) -> None:
        dialog = ContractWorkbenchDialog(self.repository_root)
        try:
            self.assertEqual(
                set(dialog.domain_buttons),
                {"formula", "geometry", "case", "runtime"},
            )
            self.assertEqual(
                dialog.domain_buttons["formula"].text(),
                "Formula & Knowledge",
            )
            self.assertEqual(dialog.domain_buttons["case"].text(), "Case")
            self.assertEqual(dialog.domain_buttons["runtime"].text(), "Runtime")
        finally:
            dialog.close()

    def test_case_domain_displays_typed_case_setup_tab(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "case")
        try:
            self.assertGreaterEqual(dialog.contract_list.count(), 1)
            self.assertEqual(dialog.tabs.tabText(dialog.tabs.count() - 1), "Case Setup")
            self.assertIn("SLAU2", dialog.case_view.toPlainText())
            self.assertIn("Blasius", dialog.case_view.toPlainText())
            case_text = dialog.case_view.toPlainText()
            self.assertIn(
                "Blasius analytical reference - skin-friction coefficient (cf)",
                case_text,
            )
            self.assertIn(
                "Blasius analytical reference - streamwise velocity profile (u)",
                case_text,
            )
            self.assertIn(
                "Blasius analytical reference - wall-normal velocity profile (v)",
                case_text,
            )
            self.assertIsNotNone(dialog.open_case_directory_button)
            self.assertIsNotNone(dialog.open_case_references_button)
            desktop_services = Mock()
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.QDesktopServices",
                desktop_services,
            ):
                dialog._open_case_directory()
            expected_directory = (
                self.repository_root / "test/cfd/plateuns2dslau2"
            ).resolve()
            desktop_services.openUrl.assert_called_once()
            self.assertEqual(
                Path(
                    desktop_services.openUrl.call_args.args[0].toLocalFile()
                ).resolve(),
                expected_directory,
            )
            desktop_services.reset_mock()
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.QDesktopServices",
                desktop_services,
            ):
                dialog._open_case_references()
            opened_paths = [
                call.args[0].toLocalFile()
                for call in desktop_services.openUrl.call_args_list
            ]
            self.assertEqual(
                [Path(path).resolve() for path in opened_paths],
                [
                    (self.repository_root / path).resolve()
                    for path in (
                        "test/suites/cpu-serial.txt",
                        "test/baselines/residual-baseline.json",
                        "test/cfd/plateuns2dslau2/results/cf.blasius",
                        "test/cfd/plateuns2dslau2/results/u.blasius",
                        "test/cfd/plateuns2dslau2/results/v.blasius",
                    )
                ],
            )
        finally:
            dialog.close()

    def test_runtime_verification_shows_resolution_separately_from_declaration(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/process-lifecycle"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            evidence_text = dialog.evidence_view.toPlainText()
            self.assertIn("[RESOLVED] google_test:", evidence_text)
            self.assertIn("Declaration: linked_not_executed", evidence_text)
            self.assertIn("Resolved file:", evidence_text)
            self.assertIn(
                "Build target: simu_imp_case_lifecycle_integration_test",
                evidence_text,
            )
            self.assertIn("Build declaration: target_resolved", evidence_text)
            self.assertIn("CTest discovery declaration: registered", evidence_text)
            self.assertIn("CTest discovery mode: PRE_TEST", evidence_text)
        finally:
            dialog.close()

    def test_gui_process_contract_displays_python_test_declarations(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/gui-process-execution"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            evidence_text = dialog.evidence_view.toPlainText()
            self.assertIn("[RESOLVED] python_test:", evidence_text)
            self.assertIn(
                "ProcessControlTest.test_terminate_and_reap_waits_for_direct_child",
                evidence_text,
            )
            self.assertIn("Declaration: linked_not_executed", evidence_text)
            self.assertIn("Resolved file:", evidence_text)
        finally:
            dialog.close()

    def test_open_declared_python_test_sources_from_verification(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            contract_index = next(
                index
                for index, summary in enumerate(dialog._summaries)
                if summary.id == "runtime/gui-process-execution"
            )
            dialog.contract_list.setCurrentRow(contract_index)
            self.assertTrue(dialog.open_tests_button.isEnabled())
            desktop_services = Mock()
            with patch(
                "tools.oneflow_gui.contract_workbench_dialog.QDesktopServices",
                desktop_services,
            ):
                dialog.open_tests_button.click()
            opened_paths = {
                Path(call.args[0].toLocalFile()).resolve()
                for call in desktop_services.openUrl.call_args_list
            }
            self.assertEqual(
                opened_paths,
                {
                    (self.repository_root / path).resolve()
                    for path in (
                        "tools/oneflow_gui/test_process_control.py",
                        "tools/oneflow_gui/test_runner_process.py",
                        "tools/oneflow_gui/test_workflow_process.py",
                    )
                },
            )
        finally:
            dialog.close()

    def test_execution_history_separates_explicit_links_from_unlinked_runs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            store_path = Path(temporary) / "execution-records.jsonl"
            store = ExecutionRecordStore(store_path)
            common = {
                "test_id": "unittest:sample",
                "framework": "unittest",
                "test_name": "SampleTest.test_contract",
                "source_path": "tools/sample_test.py",
                "command": ("python", "-m", "unittest", "tools.sample_test"),
                "outcome": "passed",
                "started_at": "2026-11-02T01:00:00+08:00",
                "finished_at": "2026-11-02T01:00:01+08:00",
                "exit_code": 0,
                "branch": "hexin",
            }
            store.append(
                ExecutionRecord(
                    **common,
                    commit="current-sha",
                    contract_ids=("runtime/gui-process-execution",),
                    artifacts=(ExecutionArtifact("log", ".oneflow/logs/run.log"),),
                )
            )
            store.append(
                ExecutionRecord(
                    **{
                        **common,
                        "test_id": "unittest:unlinked",
                        "test_name": "SampleTest.test_unlinked",
                    },
                    commit="old-sha",
                    contract_ids=(),
                )
            )
            dialog = ContractDomainDialog(
                self.repository_root,
                "runtime",
                execution_records_path=store_path,
                current_commit="current-sha",
            )
            try:
                index = next(
                    index
                    for index, summary in enumerate(dialog._summaries)
                    if summary.id == "runtime/gui-process-execution"
                )
                dialog.contract_list.setCurrentRow(index)
                text = dialog.evidence_view.toPlainText()
                linked_section = text.split("Linked executions", 1)[1].split(
                    "Recent unlinked executions", 1
                )[0]
                unlinked_section = text.split(
                    "Recent unlinked executions", 1
                )[1]
                self.assertIn("SampleTest.test_contract", linked_section)
                self.assertIn("CURRENT", linked_section)
                self.assertIn(".oneflow/logs/run.log", linked_section)
                self.assertNotIn("unlinked", linked_section)
                self.assertIn("SampleTest.test_unlinked", unlinked_section)
                self.assertIn("STALE", unlinked_section)
                self.assertIn("not attributed to this contract", text)
            finally:
                dialog.close()

    def test_runtime_domain_uses_shared_contract_navigation(self) -> None:
        dialog = ContractDomainDialog(self.repository_root, "runtime")
        try:
            self.assertGreaterEqual(dialog.contract_list.count(), 1)
            tab_names = [
                dialog.tabs.tabText(index)
                for index in range(dialog.tabs.count())
            ]
            self.assertIn("Overview", tab_names)
            self.assertIn("Principle", tab_names)
            self.assertIn("Architecture", tab_names)
            self.assertIn("Implementation / Code", tab_names)
            self.assertIn("Verification", tab_names)
            self.assertNotIn("Case Setup", tab_names)
            self.assertIsNone(dialog.open_case_directory_button)
            self.assertIsNone(dialog.open_case_references_button)
        finally:
            dialog.close()


if __name__ == "__main__":
    unittest.main()
