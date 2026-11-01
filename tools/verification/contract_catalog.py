"""Read-only compatibility views for existing verification catalogs.

This module does not own execution. It adapts existing numerical entries and
geometry JSON scenarios into a small, Qt-independent summary model.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
from typing import Any, Iterable


@dataclass(frozen=True)
class ContractEvidence:
    """A declared evidence relationship, not an execution result."""

    kind: str
    reference: str
    status: str
    scope: str


@dataclass(frozen=True)
class ContractSummary:
    """Common read-only navigation data for a verification contract."""

    id: str
    source_id: str
    title: str
    domain: str
    kind: str
    overview: str
    principle: str
    references: tuple[str, ...]
    implementation_paths: tuple[str, ...]
    documentation_paths: tuple[str, ...]
    evidence: tuple[ContractEvidence, ...]
    claims: tuple[str, ...]
    limitations: tuple[str, ...]


def _non_empty_strings(values: Iterable[Any]) -> tuple[str, ...]:
    return tuple(str(value).strip() for value in values if str(value).strip())


def _relative_path(repository_root: Path, relpath: str) -> Path:
    """Resolve a repository-relative path and reject paths escaping the root."""
    root = repository_root.resolve()
    relative = Path(relpath)
    windows_relative = PureWindowsPath(relpath)
    if (
        relative.is_absolute()
        or windows_relative.is_absolute()
        or ".." in relative.parts
        or ".." in windows_relative.parts
    ):
        raise ValueError(f"repository path must stay relative to root: {relpath}")
    candidate = (root / relative).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"repository path escapes root: {relpath}") from exc
    return candidate


def adapt_theory_entry(entry: Any) -> ContractSummary:
    """Adapt a current TheoryEntry without changing its source catalog."""
    evidence = tuple(
        ContractEvidence(
            kind="pure_api",
            reference=api,
            status="declared_not_executed",
            scope="API is named by the catalog; linked test execution is not inferred.",
        )
        for api in _non_empty_strings(entry.pure_apis)
    )
    documentation_paths = (entry.doc_relpath,) if entry.doc_relpath else ()
    if entry.csv_relpath:
        evidence += (
            ContractEvidence(
                kind="theory_data",
                reference=entry.csv_relpath,
                status="declared_not_executed",
                scope="Theory data is catalogued; comparison status is owned by its verifier.",
            ),
        )
    return ContractSummary(
        id=f"numerical/{entry.id}",
        source_id=str(entry.id),
        title=str(entry.title),
        domain="numerical",
        kind="numerical_formula",
        overview=str(entry.background),
        principle=str(entry.principle),
        references=_non_empty_strings(entry.references),
        implementation_paths=_non_empty_strings(entry.source_relpaths),
        documentation_paths=documentation_paths,
        evidence=evidence,
        claims=_non_empty_strings((entry.proves,)),
        limitations=_non_empty_strings((entry.does_not_prove,)),
    )


def adapt_geometry_scenario(
    scenario: dict[str, Any], source_path: str | None = None
) -> ContractSummary:
    """Adapt a fixture-backed geometry scenario without running its pipeline."""
    scenario_id = str(scenario.get("id") or "").strip()
    if not scenario_id:
        raise ValueError("geometry scenario requires a non-empty id")
    title = str(scenario.get("title") or scenario_id).strip()
    description = str(scenario.get("description") or "").strip()
    tests = []
    for key in ("test", "integrationTest"):
        test = scenario.get(key)
        if not isinstance(test, dict):
            continue
        suite = str(test.get("suite") or "").strip()
        name = str(test.get("name") or "").strip()
        source = str(test.get("source") or "").strip()
        if suite and name:
            tests.append(
                ContractEvidence(
                    kind="linked_test",
                    reference=f"{suite}.{name}",
                    status="linked_not_executed",
                    scope=source or "Source path not declared.",
                )
            )
    fixture_status = "fixture_available"
    if scenario.get("result") is not None:
        fixture_status = "result_recorded_unverified"
    evidence = (
        ContractEvidence(
            kind="fixture",
            reference=source_path or scenario_id,
            status=fixture_status,
            scope="Fixture data only; this adapter does not execute production C++.",
        ),
        *tests,
    )
    notes = _non_empty_strings(scenario.get("notes") or ())
    return ContractSummary(
        id=f"geometry/{scenario_id}",
        source_id=scenario_id,
        title=title,
        domain="geometry",
        kind="scenario",
        overview=description,
        principle="See scenario notes and linked production test for the intended invariant.",
        references=(),
        implementation_paths=(),
        documentation_paths=(),
        evidence=evidence,
        claims=notes,
        limitations=(
            "Fixture-backed scenario data is not evidence that production C++ executed.",
            "Linked GoogleTest declarations are references, not recorded test results.",
        ),
    )


def validate_contract_summaries(
    repository_root: Path,
    summaries: Iterable[ContractSummary],
) -> list[str]:
    """Validate stable IDs and declared repository paths for adapted summaries."""
    issues: list[str] = []
    seen_ids: set[str] = set()
    for summary in summaries:
        if not summary.id.strip():
            issues.append("contract has an empty id")
        elif summary.id in seen_ids:
            issues.append(f"duplicate contract id: {summary.id}")
        seen_ids.add(summary.id)
        declared_paths = [
            *summary.implementation_paths,
            *summary.documentation_paths,
            *(
                item.reference
                for item in summary.evidence
                if item.kind in ("theory_data", "fixture")
                and "/" in item.reference
            ),
            *(
                item.scope
                for item in summary.evidence
                if item.kind == "linked_test"
                and item.scope.endswith((".cpp", ".cc", ".cxx"))
            ),
        ]
        for relpath in declared_paths:
            try:
                path = _relative_path(repository_root, relpath)
            except ValueError as exc:
                issues.append(f"{summary.id}: {exc}")
                continue
            if not path.is_file():
                issues.append(f"{summary.id}: missing repository file: {relpath}")
    return issues

def load_contract_catalog(
    repository_root: Path,
) -> tuple[tuple[ContractSummary, ...], tuple[str, ...]]:
    """Load all current numerical and geometry entries through read-only adapters.

    Invalid source records are reported individually so one malformed scenario
    does not hide the rest of the catalog.
    """
    from .curve_catalog import list_theory_entries
    from tools.geometry_workbench.scenarios import (
        list_scenario_paths,
        load_scenario,
    )

    root = repository_root.resolve()
    summaries: list[ContractSummary] = []
    issues: list[str] = []

    for entry in list_theory_entries():
        try:
            summaries.append(adapt_theory_entry(entry))
        except (AttributeError, TypeError, ValueError) as exc:
            issues.append(f"numerical/{getattr(entry, 'id', '<unknown>')}: {exc}")

    for path in list_scenario_paths(root):
        relative_path = path.relative_to(root).as_posix()
        try:
            scenario = load_scenario(path)
            summaries.append(adapt_geometry_scenario(scenario, relative_path))
        except (OSError, ValueError, TypeError) as exc:
            issues.append(f"{relative_path}: {exc}")

    issues.extend(validate_contract_summaries(root, summaries))
    return tuple(summaries), tuple(issues)
