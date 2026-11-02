"""Capture durable evidence for completed Python unittest processes."""

from __future__ import annotations

import hashlib
import platform
import re
import subprocess
import uuid
from pathlib import Path

from .test_runner import TestCommand
from .execution_records import ExecutionArtifact, ExecutionRecord, ExecutionRecordStore


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
    value = re.sub(r"\\.git$", "", value)
    match = re.search(r"(?:github\\.com[:/])([^/]+/[^/]+)$", value)
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
    return (
        len(arguments) >= 3
        and arguments[1:3] == ("-m", "unittest")
    ) or (
        len(arguments) >= 4
        and arguments[1:3] == ("-m", "unittest")
        and "discover" in arguments[3:]
    )


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
    log_bytes = ("\\n".join(output_lines) + ("\\n" if output_lines else "")).encode(
        "utf-8"
    )
    log_path.write_bytes(log_bytes)
    try:
        artifact_path = log_path.resolve().relative_to(root).as_posix()
    except ValueError:
        artifact_path = str(log_path.resolve())

    remote = _git_value(root, "remote", "get-url", "origin")
    record = ExecutionRecord(
        test_id="unittest-suite:" + re.sub(
            r"[^a-z0-9]+", "-", command.label.lower()
        ).strip("-"),
        framework="unittest",
        test_name=command.label,
        source_path=_unittest_source_path(command),
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
        contract_ids=(),
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
