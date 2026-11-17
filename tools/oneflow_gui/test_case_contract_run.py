"""Tests for Case contract residual-runner preparation and capture linkage."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from tools.oneflow_gui.case_contract_run import (
    CaseContractRunRequest,
    build_case_runner_command,
    case_suite_entry,
    prepare_single_case_suite,
    run_case_contract_and_capture,
)
from tools.verification.execution_records import ExecutionRecordStore


class CaseContractRunHelperTest(unittest.TestCase):
    def test_case_suite_entry_strips_test_prefix(self) -> None:
        self.assertEqual(
            case_suite_entry("test/cfd/plateuns2dslau2"),
            "cfd/plateuns2dslau2",
        )
        self.assertEqual(
            case_suite_entry("cfd/plateuns2dslau2"),
            "cfd/plateuns2dslau2",
        )

    def test_prepare_single_case_suite_writes_one_entry(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            request = CaseContractRunRequest(
                contract_id="case/laminar-flat-plate-slau2",
                case_id="plateuns2dslau2",
                case_path="test/cfd/plateuns2dslau2",
                suite="test/suites/cpu-serial.txt",
                residual_baseline="test/baselines/residual-baseline.json",
                residual_absolute_tolerance=1e-15,
                comparison_outputs=("results/res.dat",),
                repository_root=root,
                oneflow_executable=root / "OneFLOW",
            )
            suite = prepare_single_case_suite(request, root / "suite-dir")
            self.assertEqual(suite.read_text(encoding="utf-8"), "cfd/plateuns2dslau2\n")

    def test_build_case_runner_command_includes_baseline_and_tolerance(self) -> None:
        command = build_case_runner_command(
            python_executable="python",
            test_script=Path("test/test.py"),
            mpi_cmd="mpirun -np 1",
            oneflow_executable="OneFLOW",
            suite_file=Path("suite.txt"),
            residual_baseline=Path("baseline.json"),
            residual_tolerance=1e-15,
        )
        self.assertEqual(
            command,
            (
                "python",
                "test/test.py",
                "mpirun -np 1",
                "OneFLOW",
                "suite.txt",
                "baseline.json",
                "1e-15",
            ),
        )

    def test_run_case_contract_and_capture_passes_explicit_contract_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "test").mkdir()
            (root / "test" / "test.py").write_text("# stub\n", encoding="utf-8")
            exe = root / "OneFLOW"
            exe.write_text("binary", encoding="utf-8")
            baseline = root / "test" / "baselines" / "residual-baseline.json"
            baseline.parent.mkdir(parents=True)
            baseline.write_text("{}", encoding="utf-8")
            records_path = root / "records.jsonl"
            request = CaseContractRunRequest(
                contract_id="case/laminar-flat-plate-slau2",
                case_id="plateuns2dslau2",
                case_path="test/cfd/plateuns2dslau2",
                suite="test/suites/cpu-serial.txt",
                residual_baseline="test/baselines/residual-baseline.json",
                residual_absolute_tolerance=1e-15,
                comparison_outputs=("results/res.dat",),
                repository_root=root,
                oneflow_executable=exe,
                python_executable="python",
            )

            def fake_subprocess_run(command, **kwargs):
                # Residual runner launch (list argv) vs git capture (git ...).
                if command and command[0] == "git":
                    values = {
                        ("remote", "get-url", "origin"): (
                            "https://github.com/eric2003/OneFLOW.git"
                        ),
                        ("rev-parse", "HEAD"): "abc",
                        ("branch", "--show-current"): "hexin",
                    }
                    value = values.get(tuple(command[1:]), "")
                    return SimpleNamespace(
                        stdout=value, returncode=0 if value else 1, stderr=""
                    )
                return SimpleNamespace(
                    returncode=0, stdout="ok\n", stderr="", args=command
                )

            with patch(
                "tools.oneflow_gui.case_contract_run.subprocess.run",
                side_effect=fake_subprocess_run,
            ), patch(
                "tools.oneflow_gui.execution_capture.subprocess.run",
                side_effect=fake_subprocess_run,
            ):
                records = run_case_contract_and_capture(
                    request, records_path=records_path
                )
            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.outcome, "passed")
            self.assertEqual(record.contract_ids, ("case/laminar-flat-plate-slau2",))
            self.assertEqual(record.framework, "OneFLOW Case Regression")
            self.assertEqual(ExecutionRecordStore(records_path).list_records(), records)


if __name__ == "__main__":
    unittest.main()
