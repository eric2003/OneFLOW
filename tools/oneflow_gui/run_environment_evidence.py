"""Persist immutable runtime environment evidence for one Workbench run."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


MANIFEST_SCHEMA_VERSION = 1


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(document, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    temporary.replace(path)


class RunEnvironmentEvidence:
    """Capture environment metadata beside a run without modifying Case inputs."""

    def __init__(
        self,
        *,
        run_id: str,
        executable: Path,
        run_directory: Path,
        command_template: Sequence[str],
        working_directory: Path,
        project_directory: Path,
        execution_mode: str,
        timeout: float = 15.0,
    ) -> None:
        self.run_id = run_id
        self.executable = executable.expanduser().resolve()
        self.run_directory = run_directory.expanduser().resolve()
        self.environment_path = self.run_directory / "environment.json"
        self.manifest_path = self.run_directory / "run-manifest.json"
        self.started_at = _utc_now()
        self._manifest: dict[str, Any] = {
            "schema_version": MANIFEST_SCHEMA_VERSION,
            "run_id": run_id,
            "status": "running",
            "started_at": self.started_at,
            "finished_at": None,
            "executable": str(self.executable),
            "command_template": list(command_template),
            "working_directory": str(working_directory.expanduser().resolve()),
            "project_directory": str(project_directory.expanduser().resolve()),
            "execution_mode": execution_mode,
            "environment_report": "environment.json",
            "environment_schema_version": None,
            "environment_query_status": "pending",
            "return_code": None,
        }
        self.run_directory.mkdir(parents=True, exist_ok=True)
        self._capture(timeout)
        _write_json(self.manifest_path, self._manifest)

    def _capture(self, timeout: float) -> None:
        report: dict[str, Any]
        try:
            completed = subprocess.run(
                [
                    str(self.executable),
                    "--environment",
                    "--format",
                    "json",
                ],
                cwd=str(self.executable.parent),
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout,
                check=False,
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"OneFLOW environment query exited with code "
                    f"{completed.returncode}: "
                    f"{(completed.stderr or completed.stdout).strip()}"
                )
            report = json.loads(completed.stdout)
            if not isinstance(report, dict) or report.get("schema_version") != 1:
                raise ValueError("Unsupported OneFLOW environment JSON schema.")
            self._manifest["environment_schema_version"] = report["schema_version"]
            self._manifest["environment_query_status"] = "captured"
        except (OSError, RuntimeError, ValueError, json.JSONDecodeError, subprocess.TimeoutExpired) as error:
            report = {
                "schema_version": 1,
                "capture_status": "failed",
                "executable": str(self.executable),
                "captured_at": _utc_now(),
                "error": str(error),
            }
            self._manifest["environment_query_status"] = "failed"
            self._manifest["environment_query_error"] = str(error)
        report.setdefault("captured_at", _utc_now())
        report.setdefault("executable", str(self.executable))
        _write_json(self.environment_path, report)

    def finish(self, *, status: str, return_code: int | None) -> None:
        """Finalize the run manifest; the environment snapshot remains immutable."""
        self._manifest["status"] = status
        self._manifest["finished_at"] = _utc_now()
        self._manifest["return_code"] = return_code
        _write_json(self.manifest_path, self._manifest)
