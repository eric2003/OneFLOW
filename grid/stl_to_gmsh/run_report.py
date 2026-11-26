"""Human-readable and machine-readable status reports for mesh-generation runs."""
from __future__ import annotations
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

STAGE_STATUSES = frozenset({"PASS", "FAIL", "NOT_RUN", "NOT_CHECKED"})
_STATUS_LABELS = {"PASS":"PASS","FAIL":"FAIL","NOT_RUN":"NOT RUN","NOT_CHECKED":"NOT CHECKED"}

def stage_result(status: str, message: str | None = None, metrics: dict[str, Any] | None = None) -> dict[str, Any]:
    """Create one validated stage result."""
    if status not in STAGE_STATUSES:
        raise ValueError(f"Unsupported run-stage status: {status!r}")
    result: dict[str, Any] = {"status": status}
    if message:
        result["message"] = message
    if metrics:
        result["metrics"] = metrics
    return result

def build_run_report(workflow: str, stages: dict[str, dict[str, Any]], exit_code: int, *,
                     metrics: dict[str, Any] | None = None,
                     artifacts: dict[str, Any] | None = None,
                     error: dict[str, str] | None = None) -> dict[str, Any]:
    """Build the stable, versioned report envelope."""
    for name, result in stages.items():
        if result.get("status") not in STAGE_STATUSES:
            raise ValueError(f"Invalid status for stage {name!r}")
    if exit_code == 0 and any(result.get("status") == "FAIL" for result in stages.values()):
        raise ValueError("A successful run report cannot contain a failed stage")
    return {
        "schema": "oneflow.run-report", "version": 1, "workflow": workflow,
        "status": "success" if exit_code == 0 else "failed",
        "exit_code": int(exit_code),
        "created_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "stages": stages, "metrics": metrics or {}, "artifacts": artifacts or {},
        "error": error,
    }

def write_run_report(path: Path, report: dict[str, Any]) -> None:
    """Atomically write JSON beside its destination."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: str | None = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp",
                delete=False) as stream:
            temporary_path = stream.name
            json.dump(report, stream, indent=2, sort_keys=True)
            stream.write("\n")
        os.replace(temporary_path, destination)
        temporary_path = None
    finally:
        if temporary_path is not None:
            try:
                os.unlink(temporary_path)
            except FileNotFoundError:
                pass

def format_run_summary(report: dict[str, Any]) -> str:
    """Render a short terminal summary; detailed data remains in JSON."""
    success = report.get("status") == "success"
    lines = [f"=== {report.get('workflow', 'OneFLOW run')}: {'SUCCESS' if success else 'FAILED'} ===",
             f"Overall status : {'PASS' if success else 'FAIL'}", "Stage results:"]
    for name, result in report.get("stages", {}).items():
        status = result.get("status", "NOT_RUN")
        message = result.get("message")
        suffix = f" — {message}" if message else ""
        lines.append(f"  [{_STATUS_LABELS.get(status, status)}] {name}{suffix}")
    metrics = report.get("metrics", {})
    if metrics:
        lines.append("Key metrics:")
        for name, value in metrics.items():
            lines.append(f"  {name}: {value}")
    artifacts = report.get("artifacts", {})
    if artifacts:
        lines.append("Output files:")
        for name, value in artifacts.items():
            if name == "files":
                continue
            lines.append(f"  {name}: {value}")
    error = report.get("error")
    if error:
        lines.append(f"Error: {error.get('message', 'Unknown error')}")
    lines.append(f"Exit code      : {report.get('exit_code', 'unknown')}")
    return "\n".join(lines)
