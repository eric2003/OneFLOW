#!/usr/bin/env python3
"""Validate incremental OneFLOW Case lifecycle contracts."""

from __future__ import annotations

from typing import Any

SCHEMA = "oneflow.case-lifecycle"
IMPLEMENTATION_STATUSES = {
    "planned", "in_progress", "implemented", "not_applicable", "blocked",
}
VERIFICATION_STATUSES = {
    "not_run", "running", "passed", "failed", "not_applicable",
}
IMPLEMENTATION_KINDS = {
    "cpp", "python", "oneflow", "gmsh", "external", "shared",
}
ACCEPTANCE_KINDS = {
    "unit", "contract", "integration", "smoke", "physical_validation", "manual",
}


def validate_lifecycle_contract(
    contract: Any, *, expected_case_id: str | None = None,
) -> list[str]:
    """Return all structural and dependency errors without requiring future stages."""
    if not isinstance(contract, dict):
        return ["Lifecycle contract root must be an object"]

    errors: list[str] = []
    if contract.get("schema") != SCHEMA:
        errors.append(f"Unsupported lifecycle schema: {contract.get('schema')!r}")
    if isinstance(contract.get("version"), bool) or contract.get("version") != 1:
        errors.append(f"Unsupported lifecycle contract version: {contract.get('version')!r}")

    case_id = contract.get("case_id")
    if not isinstance(case_id, str) or not case_id:
        errors.append("Lifecycle contract must declare a non-empty case_id")
    elif expected_case_id is not None and case_id != expected_case_id:
        errors.append(
            f"Lifecycle case_id {case_id!r} does not match catalog Case {expected_case_id!r}"
        )
    if (not isinstance(contract.get("contract_version"), int) or isinstance(contract.get("contract_version"), bool) or contract.get("contract_version", 0) < 1):
        errors.append("Lifecycle contract_version must be a positive integer")
    if not isinstance(contract.get("verification_scope"), str) or not contract.get("verification_scope"):
        errors.append("Lifecycle verification_scope must be a non-empty string")

    stages = contract.get("stages")
    if not isinstance(stages, list) or not stages:
        return [*errors, "Lifecycle stages must be a non-empty list"]

    stage_by_id: dict[str, dict[str, Any]] = {}
    for index, stage in enumerate(stages):
        prefix = f"Lifecycle stage at index {index}"
        if not isinstance(stage, dict):
            errors.append(f"{prefix} must be an object")
            continue
        stage_id = stage.get("id")
        if not isinstance(stage_id, str) or not stage_id:
            errors.append(f"{prefix} must have a non-empty id")
            continue
        if stage_id in stage_by_id:
            errors.append(f"Duplicate lifecycle stage id: {stage_id}")
            continue
        stage_by_id[stage_id] = stage
        if not isinstance(stage.get("type"), str) or not stage.get("type"):
            errors.append(f"{stage_id}: type must be a non-empty string")
        implementation_status = stage.get("implementation_status")
        if not isinstance(implementation_status, str) or implementation_status not in IMPLEMENTATION_STATUSES:
            errors.append(f"{stage_id}: invalid implementation_status {implementation_status!r}")
        verification_status = stage.get("verification_status")
        if not isinstance(verification_status, str) or verification_status not in VERIFICATION_STATUSES:
            errors.append(f"{stage_id}: invalid verification_status {verification_status!r}")

        for field in ("depends_on", "implementations", "inputs", "outputs", "acceptance"):
            if not isinstance(stage.get(field), list):
                errors.append(f"{stage_id}: {field} must be a list")

        implementations = stage.get("implementations", [])
        if isinstance(implementations, list):
            implementation_ids: set[str] = set()
            for implementation in implementations:
                if not isinstance(implementation, dict):
                    errors.append(f"{stage_id}: implementation entries must be objects")
                    continue
                implementation_id = implementation.get("id")
                if not isinstance(implementation_id, str) or not implementation_id:
                    errors.append(f"{stage_id}: implementation must have a non-empty id")
                elif implementation_id in implementation_ids:
                    errors.append(f"{stage_id}: duplicate implementation id {implementation_id!r}")
                else:
                    implementation_ids.add(implementation_id)
                if not isinstance(implementation.get("kind"), str) or implementation.get("kind") not in IMPLEMENTATION_KINDS:
                    errors.append(f"{stage_id}: invalid implementation kind {implementation.get('kind')!r}")
                if not isinstance(implementation.get("status"), str) or implementation.get("status") not in IMPLEMENTATION_STATUSES:
                    errors.append(f"{stage_id}: invalid implementation status {implementation.get('status')!r}")

        acceptance = stage.get("acceptance", [])
        if isinstance(acceptance, list):
            acceptance_ids: set[str] = set()
            for check in acceptance:
                if not isinstance(check, dict):
                    errors.append(f"{stage_id}: acceptance entries must be objects")
                    continue
                check_id = check.get("id")
                if not isinstance(check_id, str) or not check_id:
                    errors.append(f"{stage_id}: acceptance check must have a non-empty id")
                elif check_id in acceptance_ids:
                    errors.append(f"{stage_id}: duplicate acceptance id {check_id!r}")
                else:
                    acceptance_ids.add(check_id)
                if not isinstance(check.get("kind"), str) or check.get("kind") not in ACCEPTANCE_KINDS:
                    errors.append(f"{stage_id}: invalid acceptance kind {check.get('kind')!r}")
                if not isinstance(check.get("description"), str) or not check.get("description"):
                    errors.append(f"{stage_id}: acceptance check {check_id!r} needs a description")

        for field in ("inputs", "outputs"):
            artifacts = stage.get(field, [])
            if isinstance(artifacts, list):
                artifact_ids: set[str] = set()
                for artifact in artifacts:
                    if not isinstance(artifact, dict):
                        errors.append(f"{stage_id}: {field} entries must be objects")
                        continue
                    artifact_id = artifact.get("id")
                    if not isinstance(artifact_id, str) or not artifact_id:
                        errors.append(f"{stage_id}: {field} entries need a non-empty id")
                    elif artifact_id in artifact_ids:
                        errors.append(f"{stage_id}: duplicate {field} artifact id {artifact_id!r}")
                    else:
                        artifact_ids.add(artifact_id)
                    if not isinstance(artifact.get("format"), str) or not artifact.get("format"):
                        errors.append(f"{stage_id}: artifact {artifact_id!r} needs a format")

    for stage_id, stage in stage_by_id.items():
        dependencies = stage.get("depends_on", [])
        if not isinstance(dependencies, list):
            continue
        for dependency in dependencies:
            if dependency not in stage_by_id:
                errors.append(f"{stage_id}: unknown dependency {dependency!r}")
        if stage.get("verification_status") == "passed":
            if stage.get("implementation_status") != "implemented":
                errors.append(f"{stage_id}: passed verification requires implemented stage")
            if not stage.get("acceptance"):
                errors.append(f"{stage_id}: passed verification requires acceptance checks")

    # A DFS detects cycles and verifies the prerequisite state of passed stages.
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(stage_id: str) -> None:
        if stage_id in visiting:
            errors.append(f"Lifecycle dependency cycle includes {stage_id!r}")
            return
        if stage_id in visited:
            return
        visiting.add(stage_id)
        stage = stage_by_id[stage_id]
        dependencies = stage.get("depends_on", [])
        if isinstance(dependencies, list):
            for dependency in dependencies:
                if dependency not in stage_by_id:
                    continue
                visit(dependency)
                if stage.get("verification_status") == "passed":
                    dependency_status = stage_by_id[dependency].get("verification_status")
                    if dependency_status not in {"passed", "not_applicable"}:
                        errors.append(
                            f"{stage_id}: cannot be passed while dependency {dependency!r} "
                            f"is {dependency_status!r}"
                        )
        visiting.remove(stage_id)
        visited.add(stage_id)

    for stage_id in stage_by_id:
        visit(stage_id)
    return errors
