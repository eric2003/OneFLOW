"""Conservative Contract verification status from exact execution evidence.

This module is deliberately Qt-independent. It evaluates recorded executions;
it never launches tests, changes declarations, or infers unrecorded links.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .contract_catalog import ContractSummary, EvidenceResolution
from .execution_records import ExecutionRecord


_TEST_EVIDENCE_KINDS = {"linked_test", "google_test", "python_test"}
_FINAL_STATUSES = {"passed", "failed", "blocked", "stale", "not_run", "partially_verified", "unknown"}


@dataclass(frozen=True, slots=True)
class EvidenceVerification:
    kind: str
    reference: str
    status: str
    current_records: int = 0
    stale_records: int = 0


@dataclass(frozen=True, slots=True)
class ContextVerification:
    platform: str
    configuration: str
    status: str
    passed: int
    failed: int
    blocked: int
    missing: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ContractVerificationStatus:
    contract_id: str
    declared_status: str
    status: str
    evidence: tuple[EvidenceVerification, ...]
    contexts: tuple[ContextVerification, ...]
    diagnostic: str = ""


def _normal_path(value: str, root: Path) -> str:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    try:
        relative = path.resolve().relative_to(root.resolve())
    except (OSError, ValueError):
        return ""
    return relative.as_posix().casefold()


def _record_matches(
    record: ExecutionRecord,
    summary: ContractSummary,
    kind: str,
    reference: str,
    resolved_path: str,
    repository_root: Path,
) -> bool:
    if summary.id not in record.contract_ids:
        return False
    if record.framework.casefold() == "ctest":
        return kind in {"linked_test", "google_test"} and record.test_name == reference
    if record.framework.casefold() == "unittest" and kind == "python_test":
        expected = _normal_path(resolved_path, repository_root)
        actual = _normal_path(record.source_path, repository_root)
        return bool(expected and actual and expected == actual)
    return False


def evaluate_contract_verification(
    summary: ContractSummary,
    resolutions: Iterable[EvidenceResolution],
    records: Iterable[ExecutionRecord],
    *,
    repository_root: Path,
    current_commit: str | None,
) -> ContractVerificationStatus:
    """Evaluate exact linked test evidence without merging platform/configuration.

    A Python unittest record represents its declared source module invocation;
    CTest records must match the exact registered test name. Only records whose
    Contract ID is explicit are considered. Contexts are kept separate.
    """
    declared = tuple(
        item for item in summary.evidence if item.kind in _TEST_EVIDENCE_KINDS
    )
    if not declared:
        return ContractVerificationStatus(
            summary.id, summary.status, "not_run", (), (),
            "No executable linked-test evidence is declared for this contract.",
        )

    resolution_map = {
        (item.evidence.kind, item.evidence.reference): item
        for item in resolutions
    }
    all_records = tuple(records)
    per_evidence: list[EvidenceVerification] = []
    matching_by_context: dict[tuple[str, str], dict[tuple[str, str], list[ExecutionRecord]]] = {}
    stale_by_evidence: dict[tuple[str, str], int] = {}
    current_by_evidence: dict[tuple[str, str], list[ExecutionRecord]] = {}
    stale_declarations = 0

    for evidence in declared:
        key = (evidence.kind, evidence.reference)
        resolution = resolution_map.get(key)
        if resolution is not None and resolution.status == "stale":
            stale_declarations += 1
            per_evidence.append(EvidenceVerification(evidence.kind, evidence.reference, "stale"))
            continue
        resolved_path = resolution.resolved_path if resolution is not None else ""
        matches = [
            record for record in all_records
            if _record_matches(
                record, summary, evidence.kind, evidence.reference,
                resolved_path, repository_root,
            )
        ]
        current_matches = [item for item in matches if item.freshness(current_commit) == "current"]
        stale_matches = [item for item in matches if item.freshness(current_commit) == "stale"]
        unknown_matches = [item for item in matches if item.freshness(current_commit) == "unknown"]
        stale_by_evidence[key] = len(stale_matches)
        current_by_evidence[key] = current_matches

        for record in matches:
            context = (record.platform or "unknown", record.configuration or "unknown")
            matching_by_context.setdefault(context, {}).setdefault(key, []).append(record)

        outcomes = {item.outcome for item in current_matches}
        if "failed" in outcomes:
            state = "failed"
        elif "blocked" in outcomes or "cancelled" in outcomes:
            state = "blocked"
        elif "passed" in outcomes:
            state = "passed"
        elif stale_matches:
            state = "stale"
        elif unknown_matches:
            state = "unknown"
        else:
            state = "not_run"
        per_evidence.append(
            EvidenceVerification(
                evidence.kind, evidence.reference, state,
                len(current_matches), len(stale_matches),
            )
        )

    contexts: list[ContextVerification] = []
    for (platform_name, configuration), by_evidence in sorted(matching_by_context.items()):
        states: list[str] = []
        missing: list[str] = []
        for evidence in declared:
            key = (evidence.kind, evidence.reference)
            resolution = resolution_map.get(key)
            if resolution is not None and resolution.status == "stale":
                states.append("stale")
                missing.append(evidence.reference + " (stale declaration)")
                continue
            matched = by_evidence.get(key, [])
            current = [item for item in matched if item.freshness(current_commit) == "current"]
            stale = [item for item in matched if item.freshness(current_commit) == "stale"]
            if any(item.outcome == "failed" for item in current):
                states.append("failed")
            elif any(item.outcome in {"blocked", "cancelled"} for item in current):
                states.append("blocked")
            elif any(item.outcome == "passed" for item in current):
                states.append("passed")
            elif stale:
                states.append("stale")
                missing.append(evidence.reference + " (stale execution)")
            else:
                states.append("not_run")
                missing.append(evidence.reference)
        passed = states.count("passed")
        failed = states.count("failed")
        blocked = states.count("blocked")
        if failed:
            context_status = "failed"
        elif blocked:
            context_status = "blocked"
        elif passed == len(declared):
            context_status = "passed"
        elif passed:
            context_status = "partially_verified"
        elif "stale" in states:
            context_status = "stale"
        else:
            context_status = "not_run"
        contexts.append(ContextVerification(
            platform_name, configuration, context_status,
            passed, failed, blocked, tuple(missing),
        ))

    context_states = [item.status for item in contexts]
    if "failed" in context_states or any(item.status == "failed" for item in per_evidence):
        overall = "failed"
    elif "blocked" in context_states or any(item.status == "blocked" for item in per_evidence):
        overall = "blocked"
    elif stale_declarations or any(item.status == "stale" for item in per_evidence):
        overall = "stale"
    elif contexts and all(item.status == "passed" for item in contexts):
        overall = "passed"
    elif any(item.status in {"passed", "partially_verified"} for item in contexts):
        overall = "partially_verified"
    elif any(item.status == "unknown" for item in per_evidence):
        overall = "unknown"
    else:
        overall = "not_run"

    diagnostic = (
        "Observed platform/configuration contexts are reported independently; "
        "no unconfigured CI matrix is assumed."
    )
    return ContractVerificationStatus(
        summary.id, summary.status, overall, tuple(per_evidence),
        tuple(contexts), diagnostic,
    )
