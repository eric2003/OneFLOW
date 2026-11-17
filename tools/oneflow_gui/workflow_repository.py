"""Versioned persistence for workflow definitions and execution history."""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from .task_types import TaskRunStatus, normalize_run_status
from .workflow_dependencies import validate_case_dependencies


SCHEMA_VERSION = 1
_WORKFLOW_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_RUN_STATES = {status.value for status in TaskRunStatus}


@dataclass(frozen=True, slots=True)
class WorkflowDefinition:
    """A Qt-independent, versioned workflow definition."""

    workflow_id: str
    name: str
    steps: tuple[dict[str, Any], ...]
    description: str = ""
    schema_version: int = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-compatible representation."""
        data = asdict(self)
        data["steps"] = list(self.steps)
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "WorkflowDefinition":
        """Parse and validate a workflow definition."""
        if not isinstance(data, dict):
            raise ValueError("Workflow definition must be a JSON object.")
        version = data.get("schema_version")
        if version != SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported workflow schema version: {version!r}; "
                f"expected {SCHEMA_VERSION}."
            )

        workflow_id = _validate_id(data.get("workflow_id"), "workflow_id")
        name = data.get("name")
        description = data.get("description", "")
        steps = data.get("steps")
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Workflow name must be a non-empty string.")
        if not isinstance(description, str):
            raise ValueError("Workflow description must be a string.")
        if not isinstance(steps, list):
            raise ValueError("Workflow steps must be a list.")
        for index, step in enumerate(steps):
            if not isinstance(step, dict):
                raise ValueError(f"Workflow step {index} must be an object.")
            step_type = step.get("type")
            if not isinstance(step_type, str) or not step_type.strip():
                raise ValueError(
                    f"Workflow step {index} must have a non-empty type."
                )
            if step_type == "case_graph":
                if index != 0 or len(steps) != 1:
                    raise ValueError(
                        "A case_graph workflow must contain exactly one case_graph step."
                    )
                cases = step.get("cases")
                if not isinstance(cases, list):
                    raise ValueError("A case_graph step must contain a cases list.")
                validate_case_dependencies(cases)
        return cls(
            workflow_id=workflow_id,
            name=name.strip(),
            description=description,
            steps=tuple(steps),
            schema_version=SCHEMA_VERSION,
        )


@dataclass(frozen=True, slots=True)
class WorkflowRunRecord:
    """Persisted outcome of one workflow execution."""

    run_id: str
    workflow_id: str
    status: str
    started_at: str
    finished_at: str | None = None
    executable: str | None = None
    steps: tuple[dict[str, Any], ...] = ()
    artifacts: tuple[str, ...] = ()
    message: str = ""
    note: str = ""
    log_path: str | None = None
    task_type: str = "workflow"

    def to_dict(self) -> dict[str, Any]:
        """Return the JSON-compatible representation."""
        data = asdict(self)
        data["steps"] = list(self.steps)
        data["artifacts"] = list(self.artifacts)
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "WorkflowRunRecord":
        """Parse and validate a run record."""
        if not isinstance(data, dict):
            raise ValueError("Workflow run record must be a JSON object.")
        run_id = _validate_id(data.get("run_id"), "run_id")
        workflow_id = _validate_id(data.get("workflow_id"), "workflow_id")
        status = normalize_run_status(data.get("status")).value
        started_at = data.get("started_at")
        if not isinstance(started_at, str) or not started_at.strip():
            raise ValueError("Run started_at must be a non-empty string.")
        finished_at = data.get("finished_at")
        executable = data.get("executable")
        message = data.get("message", "")
        note = data.get("note", "")
        log_path = data.get("log_path")
        task_type = data.get("task_type", "workflow")
        steps = data.get("steps", [])
        artifacts = data.get("artifacts", [])
        if finished_at is not None and not isinstance(finished_at, str):
            raise ValueError("Run finished_at must be a string or null.")
        if executable is not None and not isinstance(executable, str):
            raise ValueError("Run executable must be a string or null.")
        if not isinstance(message, str):
            raise ValueError("Run message must be a string.")
        if not isinstance(note, str):
            raise ValueError("Run note must be a string.")
        if log_path is not None and not isinstance(log_path, str):
            raise ValueError("Run log_path must be a string or null.")
        if not isinstance(task_type, str) or not task_type.strip():
            raise ValueError("Run task_type must be a non-empty string.")
        if not isinstance(steps, list) or not all(
            isinstance(step, dict) for step in steps
        ):
            raise ValueError("Run steps must be a list of objects.")
        if not isinstance(artifacts, list) or not all(
            isinstance(item, str) for item in artifacts
        ):
            raise ValueError("Run artifacts must be a list of strings.")
        return cls(
            run_id=run_id,
            workflow_id=workflow_id,
            status=status,
            started_at=started_at,
            finished_at=finished_at,
            executable=executable,
            steps=tuple(steps),
            artifacts=tuple(artifacts),
            message=message,
            note=note,
            log_path=log_path,
            task_type=task_type,
        )


class WorkflowRepository:
    """Store workflow definitions separately from append-only run records.

    The supplied root is application data storage, not the source-code tree.
    Deleting a workflow never deletes its run history.
    """

    def __init__(self, root: Path) -> None:
        self._root = Path(root).expanduser().resolve()
        self._definitions = self._root / "definitions"
        self._runs = self._root / "runs"

    def list_workflows(self) -> list[WorkflowDefinition]:
        """Load all valid workflow definitions in stable name order."""
        workflows = [
            self.load_workflow(path.stem)
            for path in sorted(self._definitions.glob("*.json"))
        ] if self._definitions.is_dir() else []
        return sorted(workflows, key=lambda item: (item.name.casefold(), item.workflow_id))

    def load_workflow(self, workflow_id: str) -> WorkflowDefinition:
        """Load a workflow by its stable identifier."""
        path = self._definition_path(workflow_id)
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
        workflow = WorkflowDefinition.from_dict(data)
        if workflow.workflow_id != workflow_id:
            raise ValueError("Workflow ID does not match its filename.")
        return workflow

    def save_workflow(self, workflow: WorkflowDefinition) -> None:
        """Create or replace a definition after validating its full payload."""
        validated = WorkflowDefinition.from_dict(workflow.to_dict())
        self._write_json(self._definition_path(validated.workflow_id), validated.to_dict())

    def delete_workflow(self, workflow_id: str) -> bool:
        """Delete only the definition; retain all run records."""
        path = self._definition_path(workflow_id)
        if not path.exists():
            return False
        path.unlink()
        return True

    def save_run(self, run: WorkflowRunRecord) -> None:
        """Persist one run record independently of its workflow definition."""
        validated = WorkflowRunRecord.from_dict(run.to_dict())
        self._write_json(self._run_path(validated.run_id), validated.to_dict())

    def list_runs(
        self,
        workflow_id: str | None = None,
        *,
        status: str | None = None,
        started_before: str | None = None,
    ) -> list[WorkflowRunRecord]:
        """Return newest-first history with optional workflow, status, and age filters."""
        if workflow_id is not None:
            _validate_id(workflow_id, "workflow_id")
        if status is not None and status not in _RUN_STATES:
            raise ValueError(f"Unsupported workflow run status: {status!r}.")
        if started_before is not None and not isinstance(started_before, str):
            raise ValueError("started_before must be an ISO timestamp string.")
        records: list[WorkflowRunRecord] = []
        if not self._runs.is_dir():
            return records
        for path in self._runs.glob("*.json"):
            with path.open("r", encoding="utf-8") as stream:
                record = WorkflowRunRecord.from_dict(json.load(stream))
            if record.run_id != path.stem:
                raise ValueError("Run ID does not match its filename.")
            if workflow_id is not None and record.workflow_id != workflow_id:
                continue
            if status is not None and record.status != status:
                continue
            if started_before is not None and record.started_at >= started_before:
                continue
            records.append(record)
        return sorted(records, key=lambda item: (item.started_at, item.run_id), reverse=True)

    def update_run_note(self, run_id: str, note: str) -> WorkflowRunRecord:
        """Update the user note without changing artifacts or execution metadata."""
        if not isinstance(note, str):
            raise ValueError("Run note must be a string.")
        path = self._run_path(run_id)
        with path.open("r", encoding="utf-8") as stream:
            record = WorkflowRunRecord.from_dict(json.load(stream))
        updated = replace(record, note=note)
        self.save_run(updated)
        return updated

    def delete_run(self, run_id: str) -> bool:
        """Delete only the history record; never remove the run directory or artifacts."""
        path = self._run_path(run_id)
        if not path.exists():
            return False
        path.unlink()
        return True

    def clear_runs(
        self,
        *,
        workflow_id: str | None = None,
        statuses: set[str] | None = None,
        started_before: str | None = None,
    ) -> int:
        """Delete matching history records only and return the number removed."""
        if statuses is not None:
            invalid = statuses - _RUN_STATES
            if invalid:
                raise ValueError(f"Unsupported workflow run status: {sorted(invalid)!r}.")
        records = self.list_runs(workflow_id, started_before=started_before)
        selected = [
            record for record in records
            if statuses is None or record.status in statuses
        ]
        for record in selected:
            self.delete_run(record.run_id)
        return len(selected)

    def _definition_path(self, workflow_id: str) -> Path:
        return self._definitions / f"{_validate_id(workflow_id, 'workflow_id')}.json"

    def _run_path(self, run_id: str) -> Path:
        return self._runs / f"{_validate_id(run_id, 'run_id')}.json"

    @staticmethod
    def _write_json(path: Path, data: dict[str, Any]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(data, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        temporary.replace(path)


def _validate_id(value: Any, field_name: str) -> str:
    if not isinstance(value, str) or not _WORKFLOW_ID_PATTERN.fullmatch(value):
        raise ValueError(
            f"{field_name} must start with an ASCII letter or digit and contain "
            "only ASCII letters, digits, '.', '_' or '-'."
        )
    if value in {".", ".."}:
        raise ValueError(f"{field_name} is not a valid identifier.")
    return value
