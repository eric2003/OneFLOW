"""Read-only compatibility views for existing verification catalogs.

This module does not own execution. It adapts existing numerical entries and
geometry JSON scenarios into a small, Qt-independent summary model.
"""

from __future__ import annotations

from dataclasses import dataclass
import ast
import json
import math
import re
import tokenize
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
class EvidenceResolution:
    """Static reference resolution, kept separate from execution status."""

    evidence: ContractEvidence
    status: str
    resolved_path: str = ""
    diagnostic: str = ""


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
    output_absolute_tolerance: float
    residual_relative_tolerance: float
    residual_nan_policy: str


@dataclass(frozen=True)
class WorkflowArtifactInput:
    """A workflow step's declared artifact input binding."""

    artifact_id: str
    setting: str
    parameter_file: str


@dataclass(frozen=True)
class WorkflowArtifactOutput:
    """A workflow step's declared artifact output."""

    artifact_id: str
    path: str
    parameter_file: str = ""
    setting: str = ""


@dataclass(frozen=True)
class WorkflowStepData:
    """One case node in a dependency workflow."""

    case_id: str
    directory: str
    depends_on: tuple[str, ...]
    inputs: tuple[WorkflowArtifactInput, ...]
    outputs: tuple[WorkflowArtifactOutput, ...]


@dataclass(frozen=True)
class WorkflowContractData:
    """Typed view of a declared dependency workflow."""

    workflow_id: str
    definition_path: str
    execution_mode: str
    steps: tuple[WorkflowStepData, ...]


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
    workflow_data: WorkflowContractData | None = None
    group_id: str = ""
    group_title: str = ""
    group_order: int = 100


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
    "workflow_definition",
    "theory_data",
    "fixture",
    "linked_test",
    "google_test",
    "case_suite",
    "output_manifest",
    "residual_baseline",
    "analytical_reference",
    "python_test",
}


_EVIDENCE_PATH_KINDS = {
    "theory_data",
    "workflow_definition",
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


def _safe_relative_workflow_path(value: Any, label: str) -> str:
    """Require a relative path without parent traversal."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty relative path")
    value = value.strip()
    path = Path(value)
    windows_path = PureWindowsPath(value)
    if (
        path.is_absolute()
        or windows_path.is_absolute()
        or ".." in path.parts
        or ".." in windows_path.parts
    ):
        raise ValueError(f"{label} must stay within its declared root")
    return value


def _adapt_workflow_steps(values: Any) -> tuple[WorkflowStepData, ...]:
    """Adapt the case-graph shape used by the Workflow contract."""
    if not isinstance(values, list) or not values:
        raise ValueError("workflow.steps must be a non-empty array")
    steps = []
    for index, raw_step in enumerate(values):
        if not isinstance(raw_step, dict):
            raise ValueError(f"workflow.steps[{index}] must be an object")
        case_id = raw_step.get("id")
        if not isinstance(case_id, str) or not case_id.strip():
            raise ValueError(f"workflow.steps[{index}].id must be non-empty")
        directory = _safe_relative_workflow_path(
            raw_step.get("directory"), f"workflow.steps[{index}].directory"
        )
        dependencies = raw_step.get("depends_on", [])
        if not isinstance(dependencies, list) or any(
            not isinstance(item, str) or not item.strip() for item in dependencies
        ):
            raise ValueError(
                f"workflow.steps[{index}].depends_on must be an array of non-empty strings"
            )
        raw_inputs = raw_step.get("inputs", [])
        raw_outputs = raw_step.get("outputs", [])
        if not isinstance(raw_inputs, list) or not isinstance(raw_outputs, list):
            raise ValueError(
                f"workflow.steps[{index}].inputs and outputs must be arrays"
            )
        inputs = []
        for input_index, item in enumerate(raw_inputs):
            if not isinstance(item, dict):
                raise ValueError(
                    f"workflow.steps[{index}].inputs[{input_index}] must be an object"
                )
            artifact_id = item.get("artifact")
            setting = item.get("setting")
            if not isinstance(artifact_id, str) or not artifact_id.strip():
                raise ValueError("workflow input artifact must be non-empty")
            if not isinstance(setting, str) or not setting.strip():
                raise ValueError("workflow input setting must be non-empty")
            parameter_file = _safe_relative_workflow_path(
                item.get("parameter_file"), f"workflow input {artifact_id} parameter_file"
            )
            inputs.append(
                WorkflowArtifactInput(artifact_id.strip(), setting.strip(), parameter_file)
            )
        outputs = []
        for output_index, item in enumerate(raw_outputs):
            if not isinstance(item, dict):
                raise ValueError(
                    f"workflow.steps[{index}].outputs[{output_index}] must be an object"
                )
            artifact_id = item.get("id")
            if not isinstance(artifact_id, str) or not artifact_id.strip():
                raise ValueError("workflow output id must be non-empty")
            artifact_path = _safe_relative_workflow_path(
                item.get("path"), f"workflow output {artifact_id} path"
            )
            parameter_file = item.get("parameter_file", "")
            setting = item.get("setting", "")
            if bool(parameter_file) != bool(setting):
                raise ValueError(
                    f"workflow output {artifact_id} must declare parameter_file and setting together"
                )
            if parameter_file:
                parameter_file = _safe_relative_workflow_path(
                    parameter_file, f"workflow output {artifact_id} parameter_file"
                )
            if not isinstance(setting, str):
                raise ValueError(f"workflow output {artifact_id} setting must be a string")
            outputs.append(
                WorkflowArtifactOutput(
                    artifact_id.strip(), artifact_path, parameter_file, setting.strip()
                )
            )
        steps.append(
            WorkflowStepData(
                case_id.strip(),
                directory,
                tuple(item.strip() for item in dependencies),
                tuple(inputs),
                tuple(outputs),
            )
        )
    step_ids = [step.case_id for step in steps]
    if len(set(step_ids)) != len(step_ids):
        raise ValueError("workflow step ids must be unique")
    if any(
        dependency not in set(step_ids)
        for step in steps
        for dependency in step.depends_on
    ):
        raise ValueError("workflow dependency references an unknown step")
    producers = {}
    for step in steps:
        for output in step.outputs:
            if output.artifact_id in producers:
                raise ValueError(
                    f"workflow artifact id is declared more than once: {output.artifact_id}"
                )
            producers[output.artifact_id] = step.case_id
    for step in steps:
        for item in step.inputs:
            producer_id = producers.get(item.artifact_id)
            if producer_id is None:
                raise ValueError(
                    f"workflow input references unknown artifact: {item.artifact_id}"
                )
            if producer_id == step.case_id or producer_id not in step.depends_on:
                raise ValueError(
                    f"workflow consumer {step.case_id} must depend directly on producer {producer_id}"
                )
    return tuple(steps)


def _workflow_steps_from_definition(
    definition: dict[str, Any],
) -> tuple[WorkflowStepData, ...]:
    """Read the supported single case-graph block from a workflow definition."""
    raw_steps = definition.get("steps")
    if not isinstance(raw_steps, list):
        raise ValueError("workflow definition steps must be an array")
    graphs = [
        item for item in raw_steps
        if isinstance(item, dict) and item.get("type") == "case_graph"
    ]
    if len(graphs) != 1:
        raise ValueError("workflow definition must contain exactly one case_graph step")
    return _adapt_workflow_steps(graphs[0].get("cases"))


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

    raw_group = definition.get("group", {})
    if not isinstance(raw_group, dict):
        raise ValueError("group must be an object")
    group_id = raw_group.get("id", "")
    group_title = raw_group.get("title", "")
    group_order = raw_group.get("order", 100)
    if definition["domain"].strip() == "parameter":
        if not isinstance(group_id, str) or not group_id.strip():
            raise ValueError("parameter contract requires group.id")
        if not isinstance(group_title, str) or not group_title.strip():
            raise ValueError("parameter contract requires group.title")
    if not isinstance(group_id, str) or not isinstance(group_title, str):
        raise ValueError("group.id and group.title must be strings")
    if not isinstance(group_order, int) or isinstance(group_order, bool) or group_order < 0:
        raise ValueError("group.order must be a non-negative integer")

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
        dimension = raw_case.get("dimension")
        if not isinstance(dimension, int) or isinstance(dimension, bool) or dimension not in (1, 2, 3):
            raise ValueError("case.dimension must be 1, 2, or 3")
        max_steps = raw_case.get("max_steps")
        if not isinstance(max_steps, int) or isinstance(max_steps, bool) or max_steps <= 0:
            raise ValueError("case.max_steps must be a positive integer")
        tolerance = raw_case.get("residual_absolute_tolerance")
        if (
            not isinstance(tolerance, (int, float))
            or isinstance(tolerance, bool)
            or not math.isfinite(tolerance)
            or tolerance < 0
        ):
            raise ValueError(
                "case.residual_absolute_tolerance must be finite and non-negative"
            )
        output_tolerance = raw_case.get("output_absolute_tolerance")
        if (
            not isinstance(output_tolerance, (int, float))
            or isinstance(output_tolerance, bool)
            or not math.isfinite(output_tolerance)
            or output_tolerance < 0
        ):
            raise ValueError(
                "case.output_absolute_tolerance must be finite and non-negative"
            )
        residual_relative_tolerance = raw_case.get("residual_relative_tolerance")
        if (
            not isinstance(residual_relative_tolerance, (int, float))
            or isinstance(residual_relative_tolerance, bool)
            or not math.isfinite(residual_relative_tolerance)
            or residual_relative_tolerance < 0
        ):
            raise ValueError(
                "case.residual_relative_tolerance must be finite and non-negative"
            )
        residual_nan_policy = raw_case.get("residual_nan_policy")
        if residual_nan_policy != "reject":
            raise ValueError("case.residual_nan_policy must be reject")
        for field in ("suite", "path", "residual_baseline"):
            _safe_relative_workflow_path(raw_case[field], f"case.{field}")
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
            if len(set(case_arrays[field])) != len(case_arrays[field]):
                raise ValueError(f"case.{field} must not contain duplicate paths")
            for value in case_arrays[field]:
                _safe_relative_workflow_path(value, f"case.{field} entry")
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
            dimension=dimension,
            flow_model=raw_case["flow_model"].strip(),
            flux_scheme=raw_case["flux_scheme"].strip(),
            max_steps=max_steps,
            input_paths=case_arrays["input_paths"],
            comparison_outputs=case_arrays["comparison_outputs"],
            residual_outputs=case_arrays["residual_outputs"],
            residual_baseline=raw_case["residual_baseline"].strip(),
            analytical_references=case_arrays["analytical_references"],
            residual_absolute_tolerance=float(tolerance),
            output_absolute_tolerance=float(output_tolerance),
            residual_relative_tolerance=float(residual_relative_tolerance),
            residual_nan_policy=residual_nan_policy,
        )

    workflow_data = None
    if definition["domain"].strip() == "workflow":
        if definition["kind"].strip() != "dependency_workflow":
            raise ValueError(
                "workflow domain currently requires kind dependency_workflow"
            )
        raw_workflow = definition.get("workflow")
        if not isinstance(raw_workflow, dict):
            raise ValueError("workflow contract requires a workflow object")
        workflow_id = raw_workflow.get("workflow_id")
        if not isinstance(workflow_id, str) or not workflow_id.strip():
            raise ValueError("workflow.workflow_id must be a non-empty string")
        definition_path = _safe_relative_workflow_path(
            raw_workflow.get("definition_path"), "workflow.definition_path"
        )
        execution_mode = raw_workflow.get("execution_mode")
        if execution_mode != "dependency_graph":
            raise ValueError("workflow.execution_mode must be dependency_graph")
        workflow_data = WorkflowContractData(
            workflow_id.strip(),
            definition_path,
            execution_mode,
            _adapt_workflow_steps(raw_workflow.get("steps")),
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
        workflow_data=workflow_data,
        group_id=group_id.strip(),
        group_title=group_title.strip(),
        group_order=group_order,
    )


def _mask_cpp_non_code(source: str) -> str:
    """Replace comments and literals with spaces while preserving line breaks."""
    raw_string = re.compile(
        r'R"(?P<delimiter>[^ ()\\\t\r\n]{0,16})\(.*?\)(?P=delimiter)"',
        re.DOTALL,
    )
    non_code = re.compile(
        r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\])*"|\'(?:\\.|[^\'\\])*\'',
        re.DOTALL,
    )

    def preserve_lines(match: re.Match[str]) -> str:
        return "".join("\n" if character == "\n" else " " for character in match.group())

    source = raw_string.sub(preserve_lines, source)
    return non_code.sub(preserve_lines, source)


def _linked_test_source_path(scope: str) -> str | None:
    """Extract the leading source path from a scope with optional explanation."""
    source_path = scope.split(";", 1)[0].strip()
    if source_path.endswith((".cpp", ".cc", ".cxx")):
        return source_path
    return None


def _python_test_source_path(scope: str) -> str | None:
    """Extract a Python test source path from an evidence scope."""
    source_path = scope.split(";", 1)[0].strip()
    return source_path if source_path.endswith(".py") else None


def _read_python_source(path: Path) -> str:
    """Read Python source using PEP 263 encoding detection, including UTF-8 BOM."""
    with tokenize.open(path) as source_file:
        return source_file.read()


def _read_python_numeric_constant(
    repository_root: Path,
    relative_path: str,
    constant_name: str,
) -> float:
    """Read a numeric module constant without importing executable test code."""
    path = _relative_path(repository_root, relative_path)
    source = _read_python_source(path)
    module = ast.parse(source, filename=str(path))
    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(
            isinstance(target, ast.Name) and target.id == constant_name
            for target in node.targets
        ):
            continue
        try:
            value = ast.literal_eval(node.value)
        except (ValueError, TypeError) as exc:
            raise ValueError(
                f"{relative_path}: {constant_name} must be a numeric literal"
            ) from exc
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0
        ):
            raise ValueError(
                f"{relative_path}: {constant_name} must be finite and non-negative"
            )
        return float(value)
    raise ValueError(f"{relative_path}: missing numeric constant {constant_name}")


def _has_python_test_declaration(source: str, reference: str) -> bool:
    """Check a unittest class/method declaration without importing the module."""
    suite, separator, test_name = reference.partition(".")
    if not separator or not suite or not test_name or "." in test_name:
        return False
    try:
        module = ast.parse(source)
    except SyntaxError:
        return False
    return any(
        isinstance(node, ast.ClassDef)
        and node.name == suite
        and any(
            isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef))
            and member.name == test_name
            for member in node.body
        )
        for node in module.body
    )


def resolve_contract_evidence(
    repository_root: Path,
    summary: ContractSummary,
) -> tuple[EvidenceResolution, ...]:
    """Resolve declared artifact references without inferring execution results."""
    results: list[EvidenceResolution] = []
    for item in summary.evidence:
        path_reference = item.kind in _EVIDENCE_PATH_KINDS
        if item.kind == "fixture" and not (
            "/" in item.reference
            or "\\" in item.reference
            or bool(PureWindowsPath(item.reference).suffix)
        ):
            path_reference = False
        source_path = (
            _python_test_source_path(item.scope)
            if item.kind == "python_test"
            else _linked_test_source_path(item.scope)
            if item.kind in ("linked_test", "google_test")
            else None
        )
        source_reference = source_path is not None
        if not path_reference and not source_reference:
            results.append(EvidenceResolution(item, "declared"))
            continue

        relative = item.reference if path_reference else source_path
        try:
            path = _relative_path(repository_root, relative)
        except ValueError as exc:
            results.append(EvidenceResolution(item, "stale", diagnostic=str(exc)))
            continue
        if not path.is_file():
            results.append(
                EvidenceResolution(
                    item,
                    "stale",
                    diagnostic=f"missing repository file: {relative}",
                )
            )
            continue

        if item.kind == "python_test":
            try:
                source = _read_python_source(path)
            except (OSError, UnicodeDecodeError, SyntaxError, LookupError) as exc:
                results.append(
                    EvidenceResolution(
                        item, "stale", str(path), f"cannot read Python test source: {exc}"
                    )
                )
                continue
            if not _has_python_test_declaration(source, item.reference):
                results.append(
                    EvidenceResolution(
                        item,
                        "stale",
                        str(path),
                        f"Python unittest declaration not found: {item.reference}",
                    )
                )
                continue
            results.append(EvidenceResolution(item, "resolved", str(path)))
            continue

        if source_reference:
            suite, separator, test_name = item.reference.partition(".")
            if not separator or not suite or not test_name:
                results.append(
                    EvidenceResolution(
                        item,
                        "stale",
                        str(path),
                        f"invalid linked test identity: {item.reference}",
                    )
                )
                continue
            try:
                source = _mask_cpp_non_code(path.read_text(encoding="utf-8"))
            except (OSError, UnicodeDecodeError) as exc:
                results.append(
                    EvidenceResolution(
                        item, "stale", str(path), f"cannot read linked test source: {exc}"
                    )
                )
                continue
            declaration = re.compile(
                rf"\b(?:TEST(?:_F|_P)?|TYPED_TEST(?:_P)?)\s*"
                rf"\(\s*{re.escape(suite)}\s*,\s*{re.escape(test_name)}\s*\)"
            )
            if not declaration.search(source):
                results.append(
                    EvidenceResolution(
                        item,
                        "stale",
                        str(path),
                        f"linked GoogleTest declaration not found: {item.reference}",
                    )
                )
                continue
        results.append(EvidenceResolution(item, "resolved", str(path)))
    return tuple(results)


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
                _linked_test_source_path(item.scope)
                for item in summary.evidence
                if item.kind in ("linked_test", "google_test")
                and _linked_test_source_path(item.scope) is not None
            ),
        )))
        if summary.workflow_data is not None:
            declared_paths.append(summary.workflow_data.definition_path)
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
                try:
                    runner_tolerance = _read_python_numeric_constant(
                        repository_root,
                        "test/test.py",
                        "FLOAT_ABS_TOLERANCE",
                    )
                    if runner_tolerance != case_data.output_absolute_tolerance:
                        issues.append(
                            f"{summary.id}: output tolerance does not match "
                            "test/test.py FLOAT_ABS_TOLERANCE"
                        )
                except (OSError, SyntaxError, ValueError) as exc:
                    issues.append(
                        f"{summary.id}: cannot validate output comparison policy: {exc}"
                    )
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
                    baseline_relative_tolerance = comparison.get(
                        "relative_tolerance"
                    )
                    if baseline_relative_tolerance != case_data.residual_relative_tolerance:
                        issues.append(
                            f"{summary.id}: residual relative tolerance does not match "
                            f"{case_data.residual_baseline}"
                        )
                    baseline_nan_policy = comparison.get("nan_policy")
                    if baseline_nan_policy != case_data.residual_nan_policy:
                        issues.append(
                            f"{summary.id}: residual NaN policy does not match "
                            f"{case_data.residual_baseline}"
                        )
                except (OSError, UnicodeDecodeError, ValueError, TypeError) as exc:
                    issues.append(
                        f"{summary.id}: cannot validate residual baseline "
                        f"{case_data.residual_baseline}: {exc}"
                    )
        for item in summary.evidence:
            if item.kind != "python_test":
                continue
            source_path_reference = _python_test_source_path(item.scope)
            if source_path_reference is None:
                issues.append(
                    f"{summary.id}: Python test evidence requires a .py source path"
                )
                continue
            try:
                source_path = _relative_path(repository_root, source_path_reference)
            except ValueError as exc:
                issues.append(f"{summary.id}: {exc}")
                continue
            if not source_path.is_file():
                continue
            try:
                source = _read_python_source(source_path)
            except (OSError, UnicodeDecodeError, SyntaxError, LookupError) as exc:
                issues.append(
                    f"{summary.id}: cannot read Python test source "
                    f"{source_path_reference}: {exc}"
                )
                continue
            if not _has_python_test_declaration(source, item.reference):
                issues.append(
                    f"{summary.id}: Python unittest declaration not found: "
                    f"{item.reference} in {source_path_reference}"
                )
        if summary.workflow_data is not None:
            workflow_data = summary.workflow_data
            try:
                definition_path = _relative_path(
                    repository_root, workflow_data.definition_path
                )
            except ValueError as exc:
                issues.append(f"{summary.id}: {exc}")
                definition_path = None
            if definition_path is not None and definition_path.is_file():
                try:
                    source_definition = json.loads(
                        definition_path.read_text(encoding="utf-8")
                    )
                    if source_definition.get("workflow_id") != workflow_data.workflow_id:
                        issues.append(
                            f"{summary.id}: workflow id does not match its definition"
                        )
                    source_steps = _workflow_steps_from_definition(source_definition)
                    if source_steps != workflow_data.steps:
                        issues.append(
                            f"{summary.id}: workflow definition steps do not match typed workflow payload"
                        )
                except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
                    issues.append(
                        f"{summary.id}: cannot validate workflow definition: {exc}"
                    )

        for item in summary.evidence:
            if item.kind not in ("linked_test", "google_test"):
                continue
            source_path_reference = _linked_test_source_path(item.scope)
            if source_path_reference is None:
                continue
            try:
                source_path = _relative_path(repository_root, source_path_reference)
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
                source = _read_python_source(source_path)
            except (OSError, UnicodeDecodeError, SyntaxError, LookupError) as exc:
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
