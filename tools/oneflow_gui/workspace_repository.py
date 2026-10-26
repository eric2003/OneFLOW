"""Versioned persistence for OneFLOW workspace and case execution definitions."""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 2
_WORKSPACE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_SETTING_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
_EXECUTION_MODES = {
    "single_case",
    "sequential_cases",
    "independent_cases",
    "dependency_workflow",
}


@dataclass(frozen=True, slots=True)
class WorkspaceArtifactInput:
    """Consume one declared artifact in a case parameter file."""

    artifact: str
    setting: str
    parameter_file: str

    def to_dict(self) -> dict[str, str]:
        return {
            "artifact": self.artifact,
            "setting": self.setting,
            "parameter_file": self.parameter_file,
        }


@dataclass(frozen=True, slots=True)
class WorkspaceArtifactOutput:
    """Produce one declared artifact, optionally configuring its output path."""

    artifact_id: str
    path: str
    parameter_file: str | None = None
    setting: str | None = None

    def to_dict(self) -> dict[str, str]:
        data = {"id": self.artifact_id, "path": self.path}
        if self.parameter_file is not None:
            data["parameter_file"] = self.parameter_file
            data["setting"] = self.setting or ""
        return data


@dataclass(frozen=True, slots=True)
class WorkspaceCase:
    """A case directory registered in a workspace."""

    case_id: str
    path: str
    depends_on: tuple[str, ...] = ()
    inputs: tuple[WorkspaceArtifactInput, ...] = ()
    outputs: tuple[WorkspaceArtifactOutput, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        data["depends_on"] = list(self.depends_on)
        data["inputs"] = [item.to_dict() for item in self.inputs]
        data["outputs"] = [item.to_dict() for item in self.outputs]
        return data


@dataclass(frozen=True, slots=True)
class WorkspaceDefinition:
    """A portable workspace definition stored in a .ofw file."""

    workspace_id: str
    name: str
    execution_mode: str
    cases: tuple[WorkspaceCase, ...]
    schema_version: int = SCHEMA_VERSION
    description: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "workspace_id": self.workspace_id,
            "name": self.name,
            "description": self.description,
            "execution_mode": self.execution_mode,
            "cases": [case.to_dict() for case in self.cases],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "WorkspaceDefinition":
        if not isinstance(data, dict):
            raise ValueError("Workspace definition must be a JSON object.")
        version = data.get("schema_version")
        if version not in {1, SCHEMA_VERSION}:
            raise ValueError(
                f"Unsupported workspace schema version: {version!r}; "
                f"supported versions are 1 and {SCHEMA_VERSION}."
            )

        workspace_id = data.get("workspace_id")
        if not isinstance(workspace_id, str) or not _WORKSPACE_ID_PATTERN.fullmatch(workspace_id):
            raise ValueError("workspace_id must contain only letters, digits, dots, underscores, or hyphens.")
        name = data.get("name")
        description = data.get("description", "")
        mode = data.get("execution_mode")
        raw_cases = data.get("cases")
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Workspace name must be a non-empty string.")
        if not isinstance(description, str):
            raise ValueError("Workspace description must be a string.")
        if not isinstance(mode, str) or mode not in _EXECUTION_MODES:
            raise ValueError(f"Unsupported execution mode: {mode!r}.")
        if not isinstance(raw_cases, list) or not raw_cases:
            raise ValueError("Workspace must contain at least one case.")

        cases: list[WorkspaceCase] = []
        seen_ids: set[str] = set()
        for index, raw_case in enumerate(raw_cases):
            if not isinstance(raw_case, dict):
                raise ValueError(f"Workspace case {index} must be an object.")
            case_id = raw_case.get("case_id")
            path = raw_case.get("path")
            dependencies = raw_case.get("depends_on", [])
            if version == 1 and (raw_case.get("inputs") or raw_case.get("outputs")):
                raise ValueError("Artifact bindings require workspace schema version 2.")
            if not isinstance(case_id, str) or not _WORKSPACE_ID_PATTERN.fullmatch(case_id):
                raise ValueError(f"Workspace case {index} has an invalid case_id.")
            if case_id in seen_ids:
                raise ValueError(f"Duplicate workspace case_id: {case_id}.")
            if not isinstance(path, str) or not path.strip():
                raise ValueError(f"Workspace case {case_id} must have a non-empty path.")
            if not isinstance(dependencies, list) or any(
                not isinstance(dependency, str) for dependency in dependencies
            ):
                raise ValueError(f"Dependencies for case {case_id} must be a list of case IDs.")
            if len(set(dependencies)) != len(dependencies):
                raise ValueError(f"Case {case_id} contains duplicate dependencies.")
            inputs = _parse_artifact_inputs(raw_case.get("inputs", []), case_id)
            outputs = _parse_artifact_outputs(raw_case.get("outputs", []), case_id)
            cases.append(
                WorkspaceCase(case_id, path, tuple(dependencies), inputs, outputs)
            )
            seen_ids.add(case_id)

        case_ids = {case.case_id for case in cases}
        for case in cases:
            unknown = set(case.depends_on) - case_ids
            if unknown:
                raise ValueError(
                    f"Case {case.case_id} depends on unknown case(s): {', '.join(sorted(unknown))}."
                )
            if case.case_id in case.depends_on:
                raise ValueError(f"Case {case.case_id} cannot depend on itself.")

        if mode == "single_case" and len(cases) != 1:
            raise ValueError("Single-case workspaces must contain exactly one case.")
        if mode == "independent_cases" and any(case.depends_on for case in cases):
            raise ValueError("Independent cases cannot declare dependencies.")
        if mode == "sequential_cases" and any(case.depends_on for case in cases):
            raise ValueError("Sequential order is defined by the case list, not dependencies.")
        if mode == "dependency_workflow":
            _validate_dependency_graph(cases)
            _validate_artifact_bindings(cases)
        elif any(case.inputs or case.outputs for case in cases):
            raise ValueError("Artifact bindings are supported only in dependency workflows.")

        return cls(
            workspace_id=workspace_id,
            name=name.strip(),
            description=description,
            execution_mode=mode,
            cases=tuple(cases),
            schema_version=SCHEMA_VERSION,
        )


def _parse_artifact_inputs(raw_inputs: Any, case_id: str) -> tuple[WorkspaceArtifactInput, ...]:
    if not isinstance(raw_inputs, list):
        raise ValueError(f"Inputs for case {case_id} must be a list.")
    inputs = []
    for index, item in enumerate(raw_inputs):
        if not isinstance(item, dict):
            raise ValueError(f"Input {index} for case {case_id} must be an object.")
        artifact = _required_text(item.get("artifact"), f"input artifact for case {case_id}")
        setting = _required_setting(item.get("setting"), f"input setting for case {case_id}")
        parameter_file = _required_relative_path(
            item.get("parameter_file"), f"input parameter_file for case {case_id}"
        )
        inputs.append(WorkspaceArtifactInput(artifact, setting, parameter_file))
    return tuple(inputs)


def _parse_artifact_outputs(raw_outputs: Any, case_id: str) -> tuple[WorkspaceArtifactOutput, ...]:
    if not isinstance(raw_outputs, list):
        raise ValueError(f"Outputs for case {case_id} must be a list.")
    outputs = []
    for index, item in enumerate(raw_outputs):
        if not isinstance(item, dict):
            raise ValueError(f"Output {index} for case {case_id} must be an object.")
        artifact_id = _required_text(item.get("id"), f"output id for case {case_id}")
        path = _required_relative_path(item.get("path"), f"output path for case {case_id}")
        parameter_file = item.get("parameter_file")
        setting = item.get("setting")
        if (parameter_file is None) != (setting is None):
            raise ValueError(
                f"Output '{artifact_id}' for case {case_id} must declare both parameter_file and setting."
            )
        if parameter_file is not None:
            parameter_file = _required_relative_path(
                parameter_file, f"output parameter_file for case {case_id}"
            )
            setting = _required_setting(setting, f"output setting for case {case_id}")
        outputs.append(WorkspaceArtifactOutput(artifact_id, path, parameter_file, setting))
    return tuple(outputs)


def _required_text(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string.")
    return value.strip()


def _required_setting(value: Any, label: str) -> str:
    setting = _required_text(value, label)
    if not _SETTING_IDENTIFIER_PATTERN.fullmatch(setting):
        raise ValueError(f"{label} must be a valid parameter setting name.")
    return setting


def _required_relative_path(value: Any, label: str) -> str:
    text = _required_text(value, label)
    path = Path(text)
    if path.is_absolute() or ".." in path.parts:
        raise ValueError(f"{label} must stay within the case directory.")
    return path.as_posix()


def _validate_artifact_bindings(cases: list[WorkspaceCase]) -> None:
    """Validate artifact references with the canonical execution-plan builder."""
    from .workflow_plan import build_case_graph_plan

    graph = [
        {
            "id": case.case_id,
            "directory": case.case_id,
            "depends_on": list(case.depends_on),
            "inputs": [item.to_dict() for item in case.inputs],
            "outputs": [item.to_dict() for item in case.outputs],
        }
        for case in cases
    ]
    build_case_graph_plan(graph, Path.cwd())


def _validate_dependency_graph(cases: list[WorkspaceCase]) -> None:
    """Reject dependency cycles before a workflow can be saved or launched."""
    dependencies = {case.case_id: set(case.depends_on) for case in cases}
    completed: set[str] = set()
    while dependencies:
        ready = {
            case_id for case_id, required in dependencies.items()
            if required <= completed
        }
        if not ready:
            raise ValueError("Dependency workflow contains a cycle.")
        completed.update(ready)
        for case_id in ready:
            del dependencies[case_id]


def load_workspace(path: str | Path) -> WorkspaceDefinition:
    """Load and validate a workspace definition from a .ofw file."""
    workspace_path = Path(path).expanduser().resolve()
    data = json.loads(workspace_path.read_text(encoding="utf-8"))
    return WorkspaceDefinition.from_dict(data)


def save_workspace(path: str | Path, workspace: WorkspaceDefinition) -> Path:
    """Validate and write a workspace definition as readable JSON."""
    workspace_path = Path(path).expanduser().resolve()
    if workspace_path.suffix.lower() != ".ofw":
        raise ValueError("Workspace files must use the .ofw extension.")
    validated = WorkspaceDefinition.from_dict(workspace.to_dict())
    workspace_path.parent.mkdir(parents=True, exist_ok=True)
    workspace_path.write_text(
        json.dumps(validated.to_dict(), indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return workspace_path
