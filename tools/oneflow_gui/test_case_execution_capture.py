"""Tests for durable evidence emitted by residual case regression runs."""

from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools.oneflow_gui.execution_capture import capture_case_execution
from tools.verification.execution_records import ExecutionRecordStore


class CaseExecutionCaptureTest(unittest.TestCase):
    def _git_result(self, arguments, **_kwargs):
        values = {
            ("remote", "get-url", "origin"): "https://github.com/eric2003/OneFLOW.git",
            ("rev-parse", "HEAD"): "def456",
            ("branch", "--show-current"): "hexin",
        }
        value = values.get(tuple(arguments[1:]), "")
        return SimpleNamespace(stdout=value, returncode=0 if value else 1)

    def test_captures_launched_case_with_explicit_contract_and_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            res = root / "results" / "res.dat"
            res.parent.mkdir(parents=True)
            res.write_text("1 1e-3\n", encoding="utf-8")
            aero = root / "results" / "aero.dat"
            aero.write_text("0.0\n", encoding="utf-8")
            baseline = root / "baselines" / "residual-baseline.json"
            baseline.parent.mkdir(parents=True)
            baseline.write_text("{}", encoding="utf-8")
            log = root / "case.log"
            log.write_text("case completed\n", encoding="utf-8")
            records_path = root / "verification" / "records.jsonl"
            data = {
                "case_id": "plateuns2dslau2",
                "case_name": "plateuns2dslau2",
                "suite": "cpu-serial",
                "run_id": "case-run-7",
                "started_at": "2026-11-06T10:00:00+00:00",
                "finished_at": "2026-11-06T10:00:12+00:00",
                "executable": str(root / "build" / "Release" / "OneFLOW"),
                "command": [
                    "mpirun",
                    "-np",
                    "1",
                    str(root / "build" / "Release" / "OneFLOW"),
                    "r",
                    str(root / "cfd" / "plateuns2dslau2"),
                ],
                "return_code": 0,
                "status": "succeeded",
                "files_ok": True,
                "residual_ok": True,
                "required_outputs_ok": True,
                "working_directory": str(root / "cfd" / "plateuns2dslau2"),
                "log_path": str(log),
                "residual_baseline": str(baseline),
                "artifacts": [str(res), str(aero)],
                "source_path": "test/test.py",
            }
            with patch(
                "tools.oneflow_gui.execution_capture.subprocess.run",
                side_effect=self._git_result,
            ):
                records = capture_case_execution(
                    data,
                    repository_root=root,
                    records_path=records_path,
                    contract_ids=("case/laminar-flat-plate-slau2",),
                )
            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.outcome, "passed")
            self.assertEqual(record.framework, "OneFLOW Case Regression")
            self.assertEqual(record.test_id, "case:cpu-serial:plateuns2dslau2")
            self.assertEqual(record.invocation_id, "case-run-7")
            self.assertEqual(record.contract_ids, ("case/laminar-flat-plate-slau2",))
            self.assertEqual(record.commit, "def456")
            self.assertEqual(record.branch, "hexin")
            self.assertEqual(record.repository, "eric2003/OneFLOW")
            self.assertEqual(record.configuration, "Release")
            self.assertEqual(record.command, tuple(data["command"]))
            roles = {item.role for item in record.artifacts}
            self.assertEqual(
                roles,
                {"case-log", "case-output", "residual-baseline"},
            )
            output = next(
                item for item in record.artifacts if item.role == "case-output"
                and item.path.endswith("res.dat")
            )
            self.assertEqual(
                output.sha256,
                hashlib.sha256(res.read_bytes()).hexdigest(),
            )
            self.assertEqual(
                ExecutionRecordStore(records_path).list_records(),
                records,
            )

    def test_does_not_record_unlaunched_case(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "records.jsonl"
            data = {
                "case_id": "plateuns2dslau2",
                "status": "skipped",
                "command": ["OneFLOW", "r", "case"],
            }
            self.assertEqual(
                capture_case_execution(
                    data, repository_root=root, records_path=path
                ),
                (),
            )
            self.assertFalse(path.exists())

    def test_records_comparison_failure_even_when_process_exit_code_is_zero(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            records_path = root / "records.jsonl"
            data = {
                "case_id": "plateuns2dslau2",
                "suite": "cpu-serial",
                "run_id": "case-run-fail",
                "started_at": "2026-11-06T11:00:00+00:00",
                "finished_at": "2026-11-06T11:00:05+00:00",
                "command": ["OneFLOW", "r", "case"],
                "return_code": 0,
                "status": "succeeded",
                "files_ok": True,
                "residual_ok": False,
                "failure_reasons": ["residual history mismatch on res.dat"],
            }
            with patch(
                "tools.oneflow_gui.execution_capture.subprocess.run",
                side_effect=self._git_result,
            ):
                records = capture_case_execution(
                    data,
                    repository_root=root,
                    records_path=records_path,
                    contract_ids=("case/laminar-flat-plate-slau2",),
                )
            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.outcome, "failed")
            self.assertEqual(record.exit_code, 0)
            self.assertIn(
                "residual history mismatch on res.dat",
                record.diagnostics,
            )
            self.assertEqual(
                record.contract_ids,
                ("case/laminar-flat-plate-slau2",),
            )

    def test_contract_ids_are_never_inferred_from_case_name(self) -> None:
        """Explicit contract_ids only; matching names/paths do not create links."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            records_path = root / "records.jsonl"
            data = {
                "case_id": "laminar-flat-plate-slau2",
                "case_name": "case/laminar-flat-plate-slau2",
                "command": ["OneFLOW", "r", "test/cfd/plateuns2dslau2"],
                "return_code": 0,
                "status": "succeeded",
                "files_ok": True,
                "residual_ok": True,
                "started_at": "2026-11-06T12:00:00+00:00",
                "finished_at": "2026-11-06T12:00:01+00:00",
            }
            with patch(
                "tools.oneflow_gui.execution_capture.subprocess.run",
                side_effect=self._git_result,
            ):
                records = capture_case_execution(
                    data,
                    repository_root=root,
                    records_path=records_path,
                    # Intentionally empty: no inference allowed.
                    contract_ids=(),
                )
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].contract_ids, ())


if __name__ == "__main__":
    unittest.main()
