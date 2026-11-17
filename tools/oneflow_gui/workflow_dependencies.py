"""Dependency graph semantics for case-based workflows.

Dependencies describe required results, not merely execution order. The graph
is independent of Qt and process launching.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class CaseDependencyGraph:
    """Validated case graph with a stable dependency-respecting order."""

    cases: tuple[dict[str, Any], ...]
    ordered_ids: tuple[str, ...]
    dependencies: dict[str, tuple[str, ...]]
    dependents: dict[str, tuple[str, ...]]

    def role(self, case_id: str) -> str:
        """Describe a case's graph role; a case may have multiple roles."""
        if case_id not in self.dependencies:
            raise KeyError(f"Unknown case id: {case_id}")
        has_inputs = bool(self.dependencies[case_id])
        has_outputs = bool(self.dependents[case_id])
        if has_inputs and has_outputs:
            return "consumer_and_producer"
        if has_inputs:
            return "consumer"
        if has_outputs:
            return "producer"
        return "independent"


def validate_case_dependencies(
    cases: Sequence[Mapping[str, object]],
) -> CaseDependencyGraph:
    """Validate case IDs and dependency edges, rejecting cycles and missing IDs.

    Each case has a stable id and optional depends_on list. Input order breaks
    ties between otherwise independent cases.
    """
    if not isinstance(cases, (list, tuple)) or not cases:
        raise ValueError("A case graph must contain at least one case.")

    case_by_id: dict[str, dict[str, Any]] = {}
    declaration_order: list[str] = []
    dependencies: dict[str, tuple[str, ...]] = {}

    for index, raw_case in enumerate(cases):
        if not isinstance(raw_case, Mapping):
            raise ValueError(f"Case {index} must be an object.")
        case_id = raw_case.get("id")
        if not isinstance(case_id, str) or not case_id.strip():
            raise ValueError(f"Case {index} must have a non-empty id.")
        case_id = case_id.strip()
        if case_id in case_by_id:
            raise ValueError(f"Duplicate case id: {case_id}")
        directory = raw_case.get("directory")
        if not isinstance(directory, str) or not directory.strip():
            raise ValueError(f"Case '{case_id}' must have a non-empty directory.")
        raw_dependencies = raw_case.get("depends_on", [])
        if not isinstance(raw_dependencies, (list, tuple)):
            raise ValueError(f"Case '{case_id}' depends_on must be a list.")
        normalized_dependencies: list[str] = []
        for dependency in raw_dependencies:
            if not isinstance(dependency, str) or not dependency.strip():
                raise ValueError(
                    f"Case '{case_id}' dependencies must be non-empty case IDs."
                )
            dependency = dependency.strip()
            if dependency == case_id:
                raise ValueError(f"Case '{case_id}' cannot depend on itself.")
            if dependency in normalized_dependencies:
                raise ValueError(
                    f"Case '{case_id}' has duplicate dependency '{dependency}'."
                )
            normalized_dependencies.append(dependency)

        normalized_case = dict(raw_case)
        normalized_case["id"] = case_id
        normalized_case["directory"] = directory.strip()
        normalized_case["depends_on"] = normalized_dependencies
        case_by_id[case_id] = normalized_case
        declaration_order.append(case_id)
        dependencies[case_id] = tuple(normalized_dependencies)

    dependents_lists: dict[str, list[str]] = {
        case_id: [] for case_id in declaration_order
    }
    for case_id, required_cases in dependencies.items():
        for required_id in required_cases:
            if required_id not in case_by_id:
                raise ValueError(
                    f"Case '{case_id}' depends on unknown case '{required_id}'."
                )
            dependents_lists[required_id].append(case_id)

    pending = {case_id: len(dependencies[case_id]) for case_id in declaration_order}
    ordered: list[str] = []
    ready = [case_id for case_id in declaration_order if pending[case_id] == 0]
    while ready:
        case_id = ready.pop(0)
        ordered.append(case_id)
        for dependent_id in dependents_lists[case_id]:
            pending[dependent_id] -= 1
            if pending[dependent_id] == 0:
                ready.append(dependent_id)

    if len(ordered) != len(declaration_order):
        cyclic_ids = [case_id for case_id in declaration_order if pending[case_id] > 0]
        raise ValueError(
            "Case dependencies contain a cycle involving: " + ", ".join(cyclic_ids)
        )

    return CaseDependencyGraph(
        cases=tuple(case_by_id[case_id] for case_id in declaration_order),
        ordered_ids=tuple(ordered),
        dependencies=dependencies,
        dependents={
            case_id: tuple(values) for case_id, values in dependents_lists.items()
        },
    )
