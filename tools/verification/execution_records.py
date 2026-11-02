"""Versioned, durable records for verification executions.

Declarations and static resolution live in contract_catalog. This module records
what actually ran; it does not launch tests or infer contract relationships.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
from typing import Any
from uuid import uuid4


SCHEMA_VERSION = 1
_OUTCOMES = {"passed", "failed", "blocked", "cancelled"}


def _timestamp(value: str, field_name: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty ISO timestamp")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"{field_name} must be an ISO timestamp") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{field_name} must include a timezone")
    return parsed


def _string_tuple(values: Any, field_name: str) -> tuple[str, ...]:
    if not isinstance(values, (list, tuple)):
        raise ValueError(f"{field_name} must be an array")
    result = tuple(values)
    if any(not isinstance(value, str) or not value.strip() for value in result):
        raise ValueError(f"{field_name} entries must be non-empty strings")
    return result


@dataclass(frozen=True, slots=True)
class ExecutionArtifact:
    """A file produced by an execution, optionally identified by checksum."""

    role: str
    path: str
    sha256: str = ""

    def __post_init__(self) -> None:
        if not self.role.strip() or not self.path.strip():
            raise ValueError("artifact role and path must be non-empty")
        if self.sha256 and (len(self.sha256) != 64 or any(
            char not in "0123456789abcdefABCDEF" for char in self.sha256
        )):
            raise ValueError("artifact sha256 must contain 64 hexadecimal characters")

    def to_dict(self) -> dict[str, str]:
        return {"role": self.role, "path": self.path, "sha256": self.sha256}


@dataclass(frozen=True, slots=True)
class ExecutionRecord:
    """Immutable facts about one completed test invocation."""

    test_id: str
    framework: str
    test_name: str
    source_path: str
    command: tuple[str, ...]
    outcome: str
    started_at: str
    finished_at: str
    exit_code: int | None = None
    repository: str = ""
    commit: str = ""
    branch: str = ""
    configuration: str = ""
    platform: str = ""
    runtime_version: str = ""
    source_line: int | None = None
    working_directory: str = ""
    contract_ids: tuple[str, ...] = ()
    artifacts: tuple[ExecutionArtifact, ...] = ()
    diagnostics: tuple[str, ...] = ()
    record_id: str = ""
    schema_version: int = SCHEMA_VERSION

    def __post_init__(self) -> None:
        for name in ("test_id", "framework", "test_name", "source_path"):
            if not getattr(self, name).strip():
                raise ValueError(f"{name} must be non-empty")
        if not isinstance(self.command, tuple) or not self.command or any(
            not isinstance(part, str) or not part for part in self.command
        ):
            raise ValueError("command must be a non-empty tuple of strings")
        if self.outcome not in _OUTCOMES:
            raise ValueError(f"unsupported execution outcome: {self.outcome}")
        started = _timestamp(self.started_at, "started_at")
        finished = _timestamp(self.finished_at, "finished_at")
        if finished < started:
            raise ValueError("finished_at must not precede started_at")
        if self.exit_code is not None and (
            isinstance(self.exit_code, bool) or not isinstance(self.exit_code, int)
        ):
            raise ValueError("exit_code must be an integer or null")
        if self.source_line is not None and (
            isinstance(self.source_line, bool)
            or not isinstance(self.source_line, int)
            or self.source_line < 1
        ):
            raise ValueError("source_line must be a positive integer or null")
        if self.schema_version != SCHEMA_VERSION:
            raise ValueError(f"unsupported execution record schema: {self.schema_version}")
        object.__setattr__(self, "contract_ids", _string_tuple(self.contract_ids, "contract_ids"))
        object.__setattr__(self, "diagnostics", _string_tuple(self.diagnostics, "diagnostics"))
        if not all(isinstance(item, ExecutionArtifact) for item in self.artifacts):
            raise ValueError("artifacts must contain ExecutionArtifact values")
        if not self.record_id:
            object.__setattr__(self, "record_id", str(uuid4()))

    def to_dict(self) -> dict[str, Any]:
        """Serialize a stable schema-v1 representation."""
        return {
            "schema_version": self.schema_version,
            "record_id": self.record_id,
            "repository": self.repository,
            "commit": self.commit,
            "branch": self.branch,
            "test": {
                "id": self.test_id,
                "framework": self.framework,
                "name": self.test_name,
                "source_path": self.source_path,
                "source_line": self.source_line,
            },
            "command": list(self.command),
            "configuration": self.configuration,
            "platform": self.platform,
            "runtime_version": self.runtime_version,
            "working_directory": self.working_directory,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "exit_code": self.exit_code,
            "outcome": self.outcome,
            "contract_ids": list(self.contract_ids),
            "diagnostics": list(self.diagnostics),
            "artifacts": [item.to_dict() for item in self.artifacts],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ExecutionRecord":
        """Parse and validate a persisted schema-v1 record."""
        if not isinstance(data, dict) or data.get("schema_version") != SCHEMA_VERSION:
            raise ValueError("unsupported or invalid execution record schema")
        test = data.get("test")
        if not isinstance(test, dict):
            raise ValueError("execution record test must be an object")
        artifacts = data.get("artifacts", [])
        if not isinstance(artifacts, list):
            raise ValueError("artifacts must be an array")
        try:
            parsed_artifacts = tuple(ExecutionArtifact(**item) for item in artifacts)
            return cls(
                test_id=test["id"],
                framework=test["framework"],
                test_name=test["name"],
                source_path=test["source_path"],
                source_line=test.get("source_line"),
                command=_string_tuple(data["command"], "command"),
                outcome=data["outcome"],
                started_at=data["started_at"],
                finished_at=data["finished_at"],
                exit_code=data.get("exit_code"),
                repository=data.get("repository", ""),
                commit=data.get("commit", ""),
                branch=data.get("branch", ""),
                configuration=data.get("configuration", ""),
                platform=data.get("platform", ""),
                runtime_version=data.get("runtime_version", ""),
                working_directory=data.get("working_directory", ""),
                contract_ids=_string_tuple(data.get("contract_ids", []), "contract_ids"),
                artifacts=parsed_artifacts,
                diagnostics=_string_tuple(data.get("diagnostics", []), "diagnostics"),
                record_id=data.get("record_id", ""),
                schema_version=data["schema_version"],
            )
        except (KeyError, TypeError) as exc:
            raise ValueError(f"invalid execution record fields: {exc}") from exc

    def freshness(self, current_commit: str | None) -> str:
        """Compare evidence revision without guessing when either revision is absent."""
        if not self.commit or not current_commit:
            return "unknown"
        return "current" if self.commit == current_commit else "stale"


class ExecutionRecordStore:
    """Append-only JSON Lines store; callers choose the user-data location."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)

    def append(self, record: ExecutionRecord) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        line = json.dumps(record.to_dict(), ensure_ascii=True, separators=(",", ":"))
        with self.path.open("a", encoding="utf-8", newline="\n") as stream:
            stream.write(line + "\n")

    def list_records(self) -> tuple[ExecutionRecord, ...]:
        if not self.path.exists():
            return ()
        records = []
        with self.path.open("r", encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, start=1):
                if not line.strip():
                    continue
                try:
                    records.append(ExecutionRecord.from_dict(json.loads(line)))
                except (json.JSONDecodeError, ValueError) as exc:
                    raise ValueError(
                        f"invalid execution record at line {line_number}: {exc}"
                    ) from exc
        return tuple(records)
