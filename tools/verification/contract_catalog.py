"""Read-only compatibility views for existing verification catalogs.

This module does not own execution. It adapts existing numerical entries and
geometry JSON scenarios into a small, Qt-independent summary model.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import re
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
    architecture: str
    status: str
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
        architecture="",
        status="implemented",
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
        architecture="",
        status="fixture_only",
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


_CONTRACT_STATUSES = {
    "proposed",
    "implemented",
    "partially_verified",
    "verified",
    "blocked",
    "deprecated",
    "fixture_only",
}


def adapt_declared_contract(definition: dict[str, Any]) -> ContractSummary:
    """Adapt a versioned, domain-owned contract definition."""
    if definition.get("schema_version") != 1:
        raise ValueError("unsupported contract schema_version; expected 1")

    required = (
        "id",
        "title",
        "domain",
        "kind",
        "status",
        "overview",
        "principle",
        "architecture",
    )
    missing = [
        name for name in required
        if not isinstance(definition.get(name), str)
        or not definition[name].strip()
    ]
    if missing:
        raise ValueError(
            "contract definition is missing required fields: "
            + ", ".join(missing)
        )

    status = definition["status"].strip()
    if status not in _CONTRACT_STATUSES:
        raise ValueError(f"unsupported contract status: {status}")

    for section_name in ("implementation", "documentation", "claims"):
        section = definition.get(section_name, {})
        if not isinstance(section, dict):
            raise ValueError(f"{section_name} must be an object")

    sequence_fields = (
        ("references", definition.get("references", ())),
        ("implementation.source_paths",
         definition.get("implementation", {}).get("source_paths", ())),
        ("documentation.paths",
         definition.get("documentation", {}).get("paths", ())),
        ("claims.proves", definition.get("claims", {}).get("proves", ())),
        ("claims.does_not_prove",
         definition.get("claims", {}).get("does_not_prove", ())),
        ("verification", definition.get("verification", ())),
    )
    for field_name, value in sequence_fields:
        if not isinstance(value, list):
            raise ValueError(f"{field_name} must be an array")

    evidence_items = []
    for index, item in enumerate(definition.get("verification", ())):
        if not isinstance(item, dict):
            raise ValueError(f"verification[{index}] must be an object")
        for field in ("kind", "reference", "status", "scope"):
            if not isinstance(item.get(field), str) or not item[field].strip():
                raise ValueError(
                    f"verification[{index}] requires non-empty {field}"
                )
        evidence_items.append(
            ContractEvidence(
                kind=item["kind"].strip(),
                reference=item["reference"].strip(),
                status=item["status"].strip(),
                scope=item["scope"].strip(),
            )
        )

    contract_id = definition["id"].strip()
    if "/" not in contract_id:
        raise ValueError("declared contract id must be namespaced")

    return ContractSummary(
        id=contract_id,
        source_id=contract_id.split("/", 1)[1],
        title=definition["title"].strip(),
        domain=definition["domain"].strip(),
        kind=definition["kind"].strip(),
        overview=definition["overview"].strip(),
        principle=definition["principle"].strip(),
        architecture=definition["architecture"].strip(),
        status=status,
        references=_non_empty_strings(definition.get("references", ())),
        implementation_paths=_non_empty_strings(
            definition.get("implementation", {}).get("source_paths", ())
        ),
        documentation_paths=_non_empty_strings(
            definition.get("documentation", {}).get("paths", ())
        ),
        evidence=tuple(evidence_items),
        claims=_non_empty_strings(definition.get("claims", {}).get("proves", ())),
        limitations=_non_empty_strings(
            definition.get("claims", {}).get("does_not_prove", ())
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

        for item in summary.evidence:
            if item.kind != "linked_test" or not item.scope.endswith(
                (".cpp", ".cc", ".cxx")
            ):
                continue
            try:
                source_path = _relative_path(repository_root, item.scope)
            except ValueError:
                continue
            if not source_path.is_file():
                continue
            suite, separator, test_name = item.reference.partition(".")
            if not separator or not suite or not test_name:
                issues.append(
                    f"{summary.id}: invalid linked test identity: {item.reference}"
                )
                continue
            try:
                source = source_path.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError) as exc:
                issues.append(
                    f"{summary.id}: cannot read linked test source {item.scope}: {exc}"
                )
                continue
            declaration = re.compile(
                rf"\bTEST(?:_F|_P)?\s*\(\s*{re.escape(suite)}\s*,"
                rf"\s*{re.escape(test_name)}\s*\)"
            )
            if not declaration.search(source):
                issues.append(
                    f"{summary.id}: linked GoogleTest declaration not found: "
                    f"{item.reference} in {item.scope}"
                )
    return issues

def load_contract_catalog(
    repository_root: Path,
) -> tuple[tuple[ContractSummary, ...], tuple[str, ...]]:
    """Load all current numerical and geometry entries through read-only adapters.

    Invalid source records are reported individually so one malformed scenario
    does not hide the rest of the catalog.
    """
    from .curve_catalog import (
        list_theory_entries,
        validate_catalog_api_references,
    )
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

    declared_root = root / "tools" / "verification" / "contracts"
    if declared_root.is_dir():
        for definition_path in sorted(declared_root.rglob("*.json")):
            relative_path = definition_path.relative_to(root).as_posix()
            try:
                definition = json.loads(
                    definition_path.read_text(encoding="utf-8")
                )
                summary = adapt_declared_contract(definition)
                summaries.append(summary)
            except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
                issues.append(f"{relative_path}: {exc}")

    for path in list_scenario_paths(root):
        relative_path = path.relative_to(root).as_posix()
        try:
            scenario = load_scenario(path)
            summaries.append(adapt_geometry_scenario(scenario, relative_path))
        except (OSError, ValueError, TypeError) as exc:
            issues.append(f"{relative_path}: {exc}")

    issues.extend(validate_contract_summaries(root, summaries))
    issues.extend(validate_catalog_api_references(root))
    return tuple(summaries), tuple(issues)
