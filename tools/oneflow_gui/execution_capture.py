"""Capture durable evidence for completed Python unittest processes."""

from __future__ import annotations

import hashlib
import platform
import re
import subprocess
import uuid
from pathlib import Path

from .test_runner import TestCommand
from ..verification.execution_records import (
    ExecutionArtifact,
    ExecutionRecord,
    ExecutionRecordStore,
)


def _git_value(repository_root: Path, *arguments: str) -> str:
    try:
        result = subprocess.run(
            ["git", *arguments],
            cwd=repository_root,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return result.stdout.strip() if result.returncode == 0 else ""


def _repository_name(remote: str) -> str:
    value = remote.strip()
    if not value:
        return ""
    value = re.sub(r"\.git$", "", value)
    match = re.search(r"(?:github\.com[:/])([^/]+/[^/]+)$", value)
    return match.group(1) if match else ""


def _unittest_source_path(command: TestCommand) -> str:
    arguments = list(command.arguments[1:])
    if "-m" in arguments:
        module_index = arguments.index("-m") + 1
        if module_index < len(arguments) and arguments[module_index] == "unittest":
            unittest_arguments = arguments[module_index + 1:]
            if unittest_arguments and unittest_arguments[0] == "discover":
                if "-s" in unittest_arguments:
                    source_index = unittest_arguments.index("-s") + 1
                    if source_index < len(unittest_arguments):
                        return Path(unittest_arguments[source_index]).as_posix()
                return "."
            if unittest_arguments:
                module = unittest_arguments[0]
                return module.replace(".", "/") + ".py"
    for argument in arguments:
        if argument.endswith(".py") and not argument.startswith("-"):
            return Path(argument).as_posix()
    return "."


def is_python_unittest_command(command: TestCommand) -> bool:
    """Return true only for commands explicitly invoking unittest."""
    arguments = command.arguments
    return len(arguments) >= 3 and arguments[1:3] == ("-m", "unittest")


def _declared_contract_ids(
    repository_root: Path,
    source_path: str,
) -> tuple[str, ...]:
    """Resolve execution links only from exact catalog-declared Python test files."""
    if not source_path.endswith(".py"):
        return ()
    try:
        from tools.verification.contract_catalog import (
            contract_ids_for_python_test_source,
            load_contract_catalog,
        )

        summaries, diagnostics = load_contract_catalog(repository_root)
        if diagnostics:
            return ()
        return contract_ids_for_python_test_source(summaries, source_path)
    except (ImportError, OSError, TypeError, ValueError, AttributeError):
        return ()


def capture_python_unittest_execution(
    command: TestCommand,
    *,
    repository_root: Path,
    records_path: Path,
    artifact_directory: Path,
    started_at: str,
    finished_at: str,
    exit_code: int,
    output_lines: tuple[str, ...] | list[str],
    cancelled: bool = False,
    configuration: str = "",
) -> ExecutionRecord:
    """Persist one record after a real unittest subprocess has finished."""
    root = Path(repository_root).resolve()
    record_id = str(uuid.uuid4())
    artifact_root = Path(artifact_directory)
    artifact_root.mkdir(parents=True, exist_ok=True)
    log_path = artifact_root / (record_id + ".log")
    log_bytes = ("\n".join(output_lines) + ("\n" if output_lines else "")).encode(
        "utf-8"
    )
    log_path.write_bytes(log_bytes)
    try:
        artifact_path = log_path.resolve().relative_to(root).as_posix()
    except ValueError:
        artifact_path = str(log_path.resolve())

    remote = _git_value(root, "remote", "get-url", "origin")
    source_path = _unittest_source_path(command)
    record = ExecutionRecord(
        test_id="unittest-suite:" + re.sub(
            r"[^a-z0-9]+", "-", command.label.lower()
        ).strip("-"),
        framework="unittest",
        test_name=command.label,
        source_path=source_path,
        command=tuple(command.arguments),
        outcome="cancelled" if cancelled else ("passed" if exit_code == 0 else "failed"),
        started_at=started_at,
        finished_at=finished_at,
        exit_code=exit_code,
        repository=_repository_name(remote),
        commit=_git_value(root, "rev-parse", "HEAD"),
        branch=_git_value(root, "branch", "--show-current"),
        configuration=configuration,
        platform=platform.platform(),
        runtime_version="Python " + platform.python_version(),
        working_directory=str(command.working_directory.resolve()),
        contract_ids=_declared_contract_ids(root, source_path),
        artifacts=(
            ExecutionArtifact(
                role="stdout",
                path=artifact_path,
                sha256=hashlib.sha256(log_bytes).hexdigest(),
            ),
        ),
        diagnostics=(),
        record_id=record_id,
    )
    ExecutionRecordStore(Path(records_path)).append(record)
    return record


def _file_artifact(role: str, value: str) -> ExecutionArtifact | None:
    """Return a checksummed artifact only when the recorded file still exists."""
    path = Path(value).expanduser()
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return ExecutionArtifact(role, str(path.resolve()), digest.hexdigest())


def capture_workflow_execution(
    run_data: dict,
    *,
    repository_root: Path,
    records_path: Path,
    contract_ids: tuple[str, ...] = (),
) -> tuple[ExecutionRecord, ...]:
    """Persist evidence for each workflow step that launched a process."""
    root = Path(repository_root).resolve()
    steps = run_data.get("steps", [])
    if not isinstance(steps, list):
        raise ValueError("workflow run steps must be an array")
    executable = str(run_data.get("executable") or "")
    configuration = Path(executable).parent.name if executable else ""
    if configuration not in {"Debug", "Release", "RelWithDebInfo", "MinSizeRel"}:
        configuration = ""
    log_path = str(run_data.get("log_path") or "")
    log_artifact = _file_artifact("workflow-log", log_path) if log_path else None
    repository = _repository_name(_git_value(root, "remote", "get-url", "origin"))
    commit = _git_value(root, "rev-parse", "HEAD")
    branch = _git_value(root, "branch", "--show-current")
    records = []
    for index, step in enumerate(steps):
        if not isinstance(step, dict):
            raise ValueError(f"workflow run step {index} must be an object")
        command = step.get("command")
        if not isinstance(command, (list, tuple)) or not command:
            continue
        if any(not isinstance(part, str) or not part for part in command):
            raise ValueError(f"workflow run step {index} command is invalid")
        status = step.get("status")
        if status not in {"succeeded", "failed", "cancelled"}:
            continue
        return_code = step.get("return_code")
        if return_code is not None and (
            isinstance(return_code, bool) or not isinstance(return_code, int)
        ):
            raise ValueError(f"workflow run step {index} return_code is invalid")
        outcome = (
            "cancelled" if status == "cancelled"
            else "passed" if status == "succeeded" and return_code == 0
            else "failed"
        )
        workflow_id = str(run_data.get("workflow_id") or "workflow")
        case_id = str(step.get("case_id") or "multi-directory")
        artifacts = [log_artifact] if log_artifact is not None else []
        raw_artifacts = step.get("artifacts", [])
        if not isinstance(raw_artifacts, list):
            raise ValueError(f"workflow run step {index} artifacts must be an array")
        for artifact_path in raw_artifacts:
            if not isinstance(artifact_path, str) or not artifact_path.strip():
                raise ValueError(f"workflow run step {index} artifact path is invalid")
            artifact = _file_artifact("workflow-output", artifact_path)
            if artifact is not None:
                artifacts.append(artifact)
        message = step.get("message", "")
        diagnostics = (message,) if outcome != "passed" and message else ()
        records.append(
            ExecutionRecord(
                test_id=f"workflow:{workflow_id}:{case_id}",
                framework="OneFLOW Workflow",
                test_name=f"{workflow_id} / {case_id}",
                source_path=str(step.get("source_path") or "tools/oneflow_gui/workflow_dialog.py"),
                command=tuple(command),
                outcome=outcome,
                started_at=step.get("started_at") or run_data["started_at"],
                finished_at=step.get("finished_at") or run_data["finished_at"],
                exit_code=return_code,
                repository=repository,
                commit=commit,
                branch=branch,
                configuration=configuration,
                platform=platform.platform(),
                runtime_version="Python " + platform.python_version(),
                working_directory=str(step.get("working_directory") or ""),
                contract_ids=contract_ids,
                artifacts=tuple(artifacts),
                diagnostics=diagnostics,
                invocation_id=str(run_data.get("run_id") or ""),
            )
        )
    store = ExecutionRecordStore(Path(records_path))
    for record in records:
        store.append(record)
    return tuple(records)
