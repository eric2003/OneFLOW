"""Launch residual-runner evidence for one Case contract without name inference."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from .execution_capture import capture_case_execution
from ..verification.execution_records import ExecutionRecord


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def case_suite_entry(case_path: str) -> str:
    """Map a repository-relative case path to a residual-suite entry."""
    normalized = case_path.replace("\\", "/").strip().strip("/")
    prefix = "test/"
    if normalized.startswith(prefix):
        normalized = normalized[len(prefix) :]
    if not normalized:
        raise ValueError("case path is empty")
    return normalized


def build_case_runner_command(
    *,
    python_executable: str,
    test_script: Path,
    mpi_cmd: str,
    oneflow_executable: str,
    suite_file: Path,
    residual_baseline: Path | None,
    residual_tolerance: float | None,
) -> tuple[str, ...]:
    """Build the argv for test/test.py for a single-case residual run."""
    command: list[str] = [
        python_executable,
        test_script.as_posix(),
        mpi_cmd,
        oneflow_executable,
        str(suite_file),
    ]
    if residual_baseline is not None:
        command.append(str(residual_baseline))
        if residual_tolerance is not None:
            command.append(f"{residual_tolerance:g}")
    return tuple(command)


@dataclass(frozen=True, slots=True)
class CaseContractRunRequest:
    """Inputs required to launch one Case contract residual run."""

    contract_id: str
    case_id: str
    case_path: str
    suite: str
    residual_baseline: str
    residual_absolute_tolerance: float
    comparison_outputs: tuple[str, ...]
    repository_root: Path
    oneflow_executable: Path
    mpi_cmd: str = ""
    python_executable: str = ""
    timeout_seconds: float | None = None


def prepare_single_case_suite(request: CaseContractRunRequest, directory: Path) -> Path:
    """Write a one-line suite file under directory and return its path."""
    directory.mkdir(parents=True, exist_ok=True)
    suite_path = directory / "case-contract-suite.txt"
    entry = case_suite_entry(request.case_path)
    suite_path.write_text(entry + "\n", encoding="utf-8")
    return suite_path


def run_case_contract_and_capture(
    request: CaseContractRunRequest,
    *,
    records_path: Path,
    log_callback=None,
) -> tuple[ExecutionRecord, ...]:
    """Run the residual runner for one Case contract and append ExecutionRecords.

    Contract links come only from ``request.contract_id``. Case names and
    directory paths never create links.
    """
    root = request.repository_root.resolve()
    executable = request.oneflow_executable.expanduser().resolve()
    if not executable.is_file():
        raise FileNotFoundError(f"OneFLOW executable not found: {executable}")

    test_script = root / "test" / "test.py"
    if not test_script.is_file():
        raise FileNotFoundError(f"Residual runner not found: {test_script}")

    baseline: Path | None = None
    if request.residual_baseline.strip():
        baseline = (root / request.residual_baseline).resolve()
        try:
            baseline.relative_to(root)
        except ValueError as error:
            raise ValueError("residual baseline is outside the repository") from error
        if not baseline.is_file():
            raise FileNotFoundError(f"Residual baseline not found: {baseline}")

    python_exe = (
        request.python_executable.strip()
        or os.environ.get("PYTHON", "").strip()
        or "python"
    )
    # Empty mpi_cmd means serial launch; test.py concatenates "mpi + exe".
    mpi_cmd = request.mpi_cmd.strip()

    run_id = f"case-contract-{uuid.uuid4().hex[:12]}"
    started_at = _utc_now()
    with tempfile.TemporaryDirectory(prefix="oneflow-case-contract-") as temporary:
        suite_path = prepare_single_case_suite(request, Path(temporary))
        # Suite paths are resolved relative to the test/ directory.
        command = build_case_runner_command(
            python_executable=python_exe,
            test_script=test_script,
            mpi_cmd=mpi_cmd,
            oneflow_executable=str(executable),
            suite_file=suite_path,
            residual_baseline=baseline,
            residual_tolerance=request.residual_absolute_tolerance,
        )

        if log_callback is not None:
            log_callback("Command: " + " ".join(command))

        try:
            completed = subprocess.run(
                list(command),
                cwd=str(root / "test"),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=request.timeout_seconds,
                check=False,
            )
            return_code = completed.returncode
            stdout = completed.stdout or ""
            stderr = completed.stderr or ""
            timed_out = False
        except subprocess.TimeoutExpired as error:
            return_code = None
            stdout = (error.stdout or "") if isinstance(error.stdout, str) else ""
            stderr = (error.stderr or "") if isinstance(error.stderr, str) else ""
            timed_out = True

        finished_at = _utc_now()
        log_path = Path(temporary) / "case-contract-run.log"
        log_path.write_text(stdout + ("\n" + stderr if stderr else ""), encoding="utf-8")

        if timed_out:
            status = "timeout"
            files_ok = False
            residual_ok = False
            failure_reasons: list[str] = ["case residual runner timed out"]
        elif return_code == 0:
            status = "succeeded"
            files_ok = True
            residual_ok = True
            failure_reasons = []
        else:
            status = "failed"
            files_ok = False
            residual_ok = False
            failure_reasons = [f"residual runner exited with code {return_code}"]

        if log_callback is not None:
            log_callback(f"Status: {status} return_code={return_code}")

        # Persist log outside the temp dir so capture can hash it.
        durable_log = root / ".oneflow-runs" / "case-contract" / run_id / "case-contract-run.log"
        durable_log.parent.mkdir(parents=True, exist_ok=True)
        durable_log.write_text(log_path.read_text(encoding="utf-8"), encoding="utf-8")

        artifacts: list[str] = []
        case_workspace = root / request.case_path
        for relative in request.comparison_outputs:
            candidate = case_workspace / relative
            if candidate.is_file():
                artifacts.append(str(candidate.resolve()))

        # ExecutionRecord forbids empty command tokens; keep serial launches representable.
        recorded_command = [part for part in command if part]
        run_data = {
            "case_id": request.case_id or case_suite_entry(request.case_path).split("/")[-1],
            "suite": Path(request.suite).name if request.suite else "",
            "run_id": run_id,
            "started_at": started_at,
            "finished_at": finished_at,
            "executable": str(executable),
            "command": recorded_command,
            "return_code": return_code,
            "status": status,
            "files_ok": files_ok,
            "residual_ok": residual_ok,
            "required_outputs_ok": files_ok,
            "working_directory": str((root / "test").resolve()),
            "log_path": str(durable_log.resolve()),
            "residual_baseline": str(baseline) if baseline is not None else "",
            "artifacts": artifacts,
            "source_path": "test/test.py",
            "failure_reasons": failure_reasons,
        }
        return capture_case_execution(
            run_data,
            repository_root=root,
            records_path=records_path,
            contract_ids=(request.contract_id,),
        )


def records_path_for_repository(repository_root: Path) -> Path:
    """Default JSONL store for Case/Workbench execution evidence."""
    path = repository_root.resolve() / ".oneflow-runs" / "verification" / "execution-records.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
