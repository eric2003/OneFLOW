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
            if execution_kind == "oneflow" else None
        )
        execution_script = (root / execution["script"]).resolve() if execution_kind == "python" else None
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
            )
        )

    if errors:
        details = "\n".join(f"- {error}" for error in errors)
        raise ValueError(f"Case catalog validation failed:\n{details}")
    return cases
