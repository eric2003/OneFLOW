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
class CaseContractData:
    """Typed CFD case payload; a case is evidence for contracts, not a contract per se."""

    suite: str
    case_path: str
    solver: str
    dimension: int
    flow_model: str
    flux_scheme: str
    max_steps: int
    input_paths: tuple[str, ...]
    comparison_outputs: tuple[str, ...]
    residual_outputs: tuple[str, ...]
    residual_baseline: str
    analytical_references: tuple[str, ...]
    residual_absolute_tolerance: float


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
    case_data: CaseContractData | None = None


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


_EVIDENCE_KINDS = {
    "pure_api",
    "theory_data",
    "fixture",
    "linked_test",
    "google_test",
    "case_suite",
    "output_manifest",
    "residual_baseline",
    "analytical_reference",
}


_EVIDENCE_PATH_KINDS = {
    "theory_data",
    "fixture",
    "case_suite",
    "output_manifest",
    "residual_baseline",
    "analytical_reference",
}


_EVIDENCE_STATUSES = {
    "declared_not_executed",
    "linked_not_executed",
    "fixture_available",
    "result_recorded_unverified",
    "not_run",
    "passed",
    "failed",
    "blocked",
    "stale",
}


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
    if not isinstance(definition, dict):
        raise ValueError("contract definition must be an object")
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
        ("references", definition.get("references", [])),
        ("implementation.source_paths",
         definition.get("implementation", {}).get("source_paths", [])),
        ("documentation.paths",
         definition.get("documentation", {}).get("paths", [])),
        ("claims.proves", definition.get("claims", {}).get("proves", [])),
        ("claims.does_not_prove",
         definition.get("claims", {}).get("does_not_prove", [])),
        ("verification", definition.get("verification", [])),
    )
    for field_name, value in sequence_fields:
        if not isinstance(value, list):
            raise ValueError(f"{field_name} must be an array")

    evidence_items = []
    for index, item in enumerate(definition.get("verification", [])):
        if not isinstance(item, dict):
            raise ValueError(f"verification[{index}] must be an object")
        for field in ("kind", "reference", "status", "scope"):
            if not isinstance(item.get(field), str) or not item[field].strip():
                raise ValueError(
                    f"verification[{index}] requires non-empty {field}"
                )
        evidence_kind = item["kind"].strip()
        if evidence_kind not in _EVIDENCE_KINDS:
            raise ValueError(
                f"verification[{index}] has unsupported kind: {evidence_kind}"
            )
        evidence_status = item["status"].strip()
        if evidence_status not in _EVIDENCE_STATUSES:
            raise ValueError(
                f"verification[{index}] has unsupported status: {evidence_status}"
            )
        evidence_items.append(
            ContractEvidence(
                kind=item["kind"].strip(),
                reference=item["reference"].strip(),
                status=evidence_status,
                scope=item["scope"].strip(),
            )
        )

    contract_id = definition["id"].strip()
    domain_prefix, separator, source_id = contract_id.partition("/")
    if not separator or not source_id or domain_prefix != definition["domain"].strip():
        raise ValueError(
            "declared contract id must be namespaced by its domain"
        )

    case_data = None
    if definition["domain"].strip() == "case":
        if definition["kind"].strip() != "cfd_regression_case":
            raise ValueError("case domain currently requires kind cfd_regression_case")
        raw_case = definition.get("case")
        if not isinstance(raw_case, dict):
            raise ValueError("case contract requires a case object")
        required_case_strings = (
            "suite",
            "path",
            "solver",
            "flow_model",
            "flux_scheme",
            "residual_baseline",
        )
        for field in required_case_strings:
            if not isinstance(raw_case.get(field), str) or not raw_case[field].strip():
                raise ValueError(f"case.{field} must be a non-empty string")
        for field in ("dimension", "max_steps"):
            if not isinstance(raw_case.get(field), int) or isinstance(raw_case[field], bool) or raw_case[field] <= 0:
                raise ValueError(f"case.{field} must be a positive integer")
        tolerance = raw_case.get("residual_absolute_tolerance")
        if not isinstance(tolerance, (int, float)) or isinstance(tolerance, bool) or tolerance < 0:
            raise ValueError("case.residual_absolute_tolerance must be non-negative")
        case_arrays = {}
        for field in (
            "input_paths",
            "comparison_outputs",
            "residual_outputs",
            "analytical_references",
        ):
            values = raw_case.get(field, [])
            if not isinstance(values, list) or any(
                not isinstance(value, str) or not value.strip() for value in values
            ):
                raise ValueError(f"case.{field} must be an array of non-empty strings")
            case_arrays[field] = _non_empty_strings(values)
        if (
            not case_arrays["input_paths"]
            or not case_arrays["comparison_outputs"]
            or not case_arrays["residual_outputs"]
        ):
            raise ValueError(
                "case requires input_paths, comparison_outputs, and residual_outputs"
            )
        if not set(case_arrays["residual_outputs"]).issubset(
            case_arrays["comparison_outputs"]
        ):
            raise ValueError(
                "case.residual_outputs must be included in comparison_outputs"
            )
        case_data = CaseContractData(
            suite=raw_case["suite"].strip(),
            case_path=raw_case["path"].strip(),
            solver=raw_case["solver"].strip(),
            dimension=raw_case["dimension"],
            flow_model=raw_case["flow_model"].strip(),
            flux_scheme=raw_case["flux_scheme"].strip(),
            max_steps=raw_case["max_steps"],
            input_paths=case_arrays["input_paths"],
            comparison_outputs=case_arrays["comparison_outputs"],
            residual_outputs=case_arrays["residual_outputs"],
            residual_baseline=raw_case["residual_baseline"].strip(),
            analytical_references=case_arrays["analytical_references"],
            residual_absolute_tolerance=float(tolerance),
        )

    return ContractSummary(
        id=contract_id,
        source_id=source_id,
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
        case_data=case_data,
    )


def _mask_cpp_non_code(source: str) -> str:
    """Replace comments and literals with spaces while preserving line breaks."""
    raw_string = re.compile(
        r'R"(?P<delimiter>[^ ()\\\\\\t\\r\\n]{0,16})\\(.*?\\)(?P=delimiter)"',
        re.DOTALL,
    )
    non_code = re.compile(
        r'//[^\\n]*|/\\*.*?\\*/|"(?:\\\\.|[^"\\\\])*"|\\'(?:\\\\.|[^\\'\\\\])*\\'',
        re.DOTALL,
    )

    def preserve_lines(match: re.Match[str]) -> str:
        return "".join("\\n" if character == "\\n" else " " for character in match.group())

    source = raw_string.sub(preserve_lines, source)
    return non_code.sub(preserve_lines, source)


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
        declared_paths = list(dict.fromkeys((
            *summary.implementation_paths,
            *summary.documentation_paths,
            *(
                item.reference
                for item in summary.evidence
                if item.kind in _EVIDENCE_PATH_KINDS
                and (
                    item.kind != "fixture"
                    or "/" in item.reference
                    or "\\" in item.reference
                    or bool(PureWindowsPath(item.reference).suffix)
                )
            ),
            *(
                item.scope
                for item in summary.evidence
                if item.kind in ("linked_test", "google_test")
                and item.scope.endswith((".cpp", ".cc", ".cxx"))
            ),
        )))
        if summary.case_data is not None:
            case_data = summary.case_data
            declared_paths.extend(
                (
                    case_data.suite,
                    *case_data.input_paths,
                    *case_data.analytical_references,
                    case_data.residual_baseline,
                )
            )
        for relpath in declared_paths:
            try:
                path = _relative_path(repository_root, relpath)
            except ValueError as exc:
                issues.append(f"{summary.id}: {exc}")
                continue
            if not path.is_file():
                issues.append(f"{summary.id}: missing repository file: {relpath}")

        if summary.case_data is not None:
            case_data = summary.case_data
            try:
                suite_path = _relative_path(repository_root, case_data.suite)
                case_path = _relative_path(repository_root, case_data.case_path)
            except ValueError:
                suite_path = None
                case_path = None
            if suite_path is not None and suite_path.is_file():
                try:
                    suite_text = suite_path.read_text(encoding="utf-8-sig")
                    expected_entry = case_data.case_path.removeprefix("test/")
                    if expected_entry not in suite_text:
                        issues.append(
                            f"{summary.id}: case path is not declared by suite "
                            f"{case_data.suite}: {expected_entry}"
                        )
                except (OSError, UnicodeDecodeError) as exc:
                    issues.append(
                        f"{summary.id}: cannot read suite {case_data.suite}: {exc}"
                    )
            if case_path is not None and not case_path.is_dir():
                issues.append(
                    f"{summary.id}: missing case directory: {case_data.case_path}"
                )
            if case_path is not None and case_path.is_dir():
                compare_path = case_path / "autotest" / "compare.txt"
                if not compare_path.is_file():
                    issues.append(
                        f"{summary.id}: missing autotest/compare.txt"
                    )
                else:
                    try:
                        declared_outputs = {
                            line.strip()
                            for line in compare_path.read_text(encoding="utf-8-sig").splitlines()
                            if line.strip() and not line.lstrip().startswith("#")
                        }
                        contract_outputs = set(case_data.comparison_outputs)
                        missing_outputs = contract_outputs - declared_outputs
                        untracked_outputs = declared_outputs - contract_outputs
                        if missing_outputs:
                            issues.append(
                                f"{summary.id}: comparison outputs are not declared "
                                f"by autotest/compare.txt: {', '.join(sorted(missing_outputs))}"
                            )
                        if untracked_outputs:
                            issues.append(
                                f"{summary.id}: autotest/compare.txt contains outputs "
                                f"missing from the contract: {', '.join(sorted(untracked_outputs))}"
                            )
                    except (OSError, UnicodeDecodeError) as exc:
                        issues.append(
                            f"{summary.id}: cannot read autotest/compare.txt: {exc}"
                        )

            try:
                residual_path = _relative_path(
                    repository_root, case_data.residual_baseline
                )
            except ValueError:
                residual_path = None
            if residual_path is not None and residual_path.is_file():
                try:
                    residual_data = json.loads(
                        residual_path.read_text(encoding="utf-8")
                    )
                    if not isinstance(residual_data, dict):
                        raise ValueError("residual baseline root must be an object")
                    if residual_data.get("schema") != "oneflow.residual-baseline":
                        issues.append(
                            f"{summary.id}: unsupported residual baseline schema"
                        )
                    baseline_cases = residual_data.get("cases", {})
                    if not isinstance(baseline_cases, dict):
                        raise ValueError("residual baseline cases must be an object")
                    case_name = case_data.case_path.removeprefix("test/").split("/")[-1]
                    baseline_case = baseline_cases.get(case_name)
                    if not isinstance(baseline_case, dict):
                        issues.append(
                            f"{summary.id}: case is missing from residual baseline "
                            f"{case_data.residual_baseline}"
                        )
                    else:
                        baseline_files = baseline_case.get("files", {})
                        for residual_output in case_data.residual_outputs:
                            baseline_file = baseline_files.get(residual_output)
                            if not isinstance(baseline_file, dict):
                                issues.append(
                                    f"{summary.id}: residual output is missing from "
                                    f"baseline: {residual_output}"
                                )
                                continue
                            variables = baseline_file.get("variables", [])
                            if not {"iter", "sub-iter"}.issubset(variables):
                                issues.append(
                                    f"{summary.id}: residual baseline lacks iteration "
                                    f"columns for {residual_output}"
                                )
                            if baseline_file.get("summary", {}).get("row_count") != case_data.max_steps:
                                issues.append(
                                    f"{summary.id}: residual row count does not match "
                                    f"case max_steps for {residual_output}"
                                )
                    comparison = residual_data.get("comparison", {})
                    if not isinstance(comparison, dict):
                        raise ValueError("residual baseline comparison must be an object")
                    baseline_tolerance = comparison.get("absolute_tolerance")
                    if baseline_tolerance != case_data.residual_absolute_tolerance:
                        issues.append(
                            f"{summary.id}: residual tolerance does not match "
                            f"{case_data.residual_baseline}"
                        )
                except (OSError, UnicodeDecodeError, ValueError, TypeError) as exc:
                    issues.append(
                        f"{summary.id}: cannot validate residual baseline "
                        f"{case_data.residual_baseline}: {exc}"
                    )
        for item in summary.evidence:
            if item.kind not in ("linked_test", "google_test") or not item.scope.endswith(
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
            source = _mask_cpp_non_code(source)
            declaration = re.compile(
                rf"\b(?:TEST(?:_F|_P)?|TYPED_TEST(?:_P)?)\s*"
                rf"\(\s*{re.escape(suite)}\s*,\s*{re.escape(test_name)}\s*\)"
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
    """Load numerical, geometry, and declared cross-domain contracts.

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
