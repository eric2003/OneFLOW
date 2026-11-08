"""Load and validate case definitions from the repository catalog."""

from __future__ import annotations

from pathlib import Path

from tools.mesh_cases import load_catalog, validate_case

from .models import CaseDefinition


def load_cases(repository_root: Path) -> list[CaseDefinition]:
    """Resolve registered cases and reject invalid catalog/profile entries."""

    root = repository_root.expanduser().resolve()
    catalog_path = root / "grid" / "case_catalog.json"
    catalog = load_catalog(catalog_path)
    cases: list[CaseDefinition] = []
    errors: list[str] = []

    for entry in catalog["cases"]:
        errors.extend(validate_case(root, entry))
        case_directory = (root / entry["directory"]).resolve()
        profile_value = entry.get("profile")
        layout_value = entry.get("layout")
        execution = entry.get("execution") or {}
        execution_kind = str(execution.get("kind") or "oneflow")
        runtime_config = (
            (case_directory / entry["runtimeConfig"]).resolve()
            if execution_kind == "oneflow" and entry.get("runtimeConfig") else None
        )
        execution_script = (root / execution["script"]).resolve() if execution_kind == "python" else None
        directory_parts = Path(entry["directory"]).parts
        workspace_relative_directory = (
            Path(*directory_parts[1:])
            if directory_parts and directory_parts[0] == "grid"
            else None
        )
        cases.append(
            CaseDefinition(
                case_id=entry["id"],
                title=entry.get("title", entry["id"]),
                case_directory=case_directory,
                runtime_config=runtime_config if execution_kind == "oneflow" else None,
                profile=(case_directory / profile_value).resolve() if profile_value else None,
                layout=(case_directory / layout_value).resolve() if layout_value else None,
                group=str(entry.get("group") or "General"),
                execution_kind=execution_kind,
                execution_script=execution_script,
                workspace_relative_directory=workspace_relative_directory,
            )
        )

    if errors:
        details = "\n".join(f"- {error}" for error in errors)
        raise ValueError(f"Case catalog validation failed:\n{details}")
    return cases


def select_cases_by_execution_kind(
    cases: list[CaseDefinition], execution_kind: str,
) -> list[CaseDefinition]:
    """Return only cases belonging to one Development Tests execution mode."""
    if execution_kind not in {"oneflow", "python"}:
        raise ValueError(f"Unsupported development case execution kind: {execution_kind!r}")
    return [case for case in cases if case.execution_kind == execution_kind]
