"""Build explicit execution plans for dependent OneFLOW cases.

The plan records artifact locations and consumer overrides without launching
OneFLOW or mutating case parameter files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from .workflow_dependencies import validate_case_dependencies


@dataclass(frozen=True, slots=True)
class ArtifactBinding:
    """Connect a producer output to a consumer's runtime input setting."""

    artifact_id: str
    producer_id: str
    producer_directory: Path
    relative_path: Path
    consumer_id: str
    consumer_setting: str
    parameter_file: Path

    @property
    def source_case_directory(self) -> Path:
        return self.producer_directory

    @property
    def source_path(self) -> Path:
        return self.producer_directory / self.relative_path

    @property
    def runtime_overrides(self) -> dict[str, str]:
        """Values needed to read this artifact from another case directory."""
        return {
            "sourceGridCaseDir": str(self.source_case_directory),
            self.consumer_setting: self.relative_path.as_posix(),
        }


@dataclass(frozen=True, slots=True)
class WorkflowPlanStep:
    """One case in dependency-respecting execution order."""

    case_id: str
    directory: Path
    depends_on: tuple[str, ...]
    inputs: tuple[ArtifactBinding, ...]


@dataclass(frozen=True, slots=True)
class WorkflowExecutionPlan:
    """Validated case order plus explicit producer-consumer artifact bindings."""

    steps: tuple[WorkflowPlanStep, ...]


def build_case_graph_plan(
    cases: Sequence[Mapping[str, Any]],
    base_directory: Path | str,
) -> WorkflowExecutionPlan:
    """Build a plan from case graph metadata.

    Each producer declares outputs as objects with id and path. Each consumer
    declares inputs as objects with artifact and setting. Paths are relative
    to the producer case directory, so grid/mesh.ofl and results/mesh.ofl
    are equally valid.
    """
    graph = validate_case_dependencies(cases)
    base = Path(base_directory).expanduser().resolve()
    case_by_id = {case["id"]: case for case in graph.cases}
    directories = {
        case_id: _resolve_case_directory(base, case["directory"], case_id)
        for case_id, case in case_by_id.items()
    }

    outputs: dict[str, tuple[str, Path]] = {}
    for case_id, case in case_by_id.items():
        raw_outputs = case.get("outputs", [])
        if not isinstance(raw_outputs, list):
            raise ValueError(f"Case '{case_id}' outputs must be a list.")
        for output in raw_outputs:
            if not isinstance(output, Mapping):
                raise ValueError(f"Case '{case_id}' output must be an object.")
            artifact_id = _non_empty_string(output.get("id"), "output id", case_id)
            relative_path = _relative_artifact_path(
                output.get("path"), f"output '{artifact_id}' in case '{case_id}'"
            )
            if artifact_id in outputs:
                raise ValueError(f"Duplicate artifact id: {artifact_id}")
            outputs[artifact_id] = (case_id, relative_path)

    inputs_by_case: dict[str, list[ArtifactBinding]] = {
        case_id: [] for case_id in case_by_id
    }
    for case_id, case in case_by_id.items():
        raw_inputs = case.get("inputs", [])
        if not isinstance(raw_inputs, list):
            raise ValueError(f"Case '{case_id}' inputs must be a list.")
        for item in raw_inputs:
            if not isinstance(item, Mapping):
                raise ValueError(f"Case '{case_id}' input must be an object.")
            artifact_id = _non_empty_string(
                item.get("artifact"), "input artifact", case_id
            )
            setting = _non_empty_string(item.get("setting"), "input setting", case_id)
            parameter_file = _relative_artifact_path(
                item.get("parameter_file"),
                f"parameter file for input '{artifact_id}' in case '{case_id}'",
            )
            if artifact_id not in outputs:
                raise ValueError(
                    f"Case '{case_id}' consumes unknown artifact '{artifact_id}'."
                )
            producer_id, relative_path = outputs[artifact_id]
            if producer_id not in graph.dependencies[case_id]:
                raise ValueError(
                    f"Case '{case_id}' must depend directly on producer "
                    f"'{producer_id}' for artifact '{artifact_id}'."
                )
            inputs_by_case[case_id].append(
                ArtifactBinding(
                    artifact_id=artifact_id,
                    producer_id=producer_id,
                    producer_directory=directories[producer_id],
                    relative_path=relative_path,
                    consumer_id=case_id,
                    consumer_setting=setting,
                    parameter_file=parameter_file,
                )
            )

    steps = tuple(
        WorkflowPlanStep(
            case_id=case_id,
            directory=directories[case_id],
            depends_on=graph.dependencies[case_id],
            inputs=tuple(inputs_by_case[case_id]),
        )
        for case_id in graph.ordered_ids
    )
    return WorkflowExecutionPlan(steps=steps)


def _resolve_case_directory(base: Path, value: str, case_id: str) -> Path:
    path = Path(value).expanduser()
    resolved = path.resolve() if path.is_absolute() else (base / path).resolve()
    if not path.is_absolute():
        try:
            resolved.relative_to(base)
        except ValueError as error:
            raise ValueError(
                f"Case '{case_id}' directory must stay within the workflow base."
            ) from error
    return resolved


def _relative_artifact_path(value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} path must be a non-empty relative path.")
    path = Path(value.strip())
    if path.is_absolute() or any(part == ".." for part in path.parts):
        raise ValueError(f"{label} path must stay within its producer case.")
    return path


def _non_empty_string(value: Any, field: str, case_id: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Case '{case_id}' {field} must be a non-empty string.")
    return value.strip()
