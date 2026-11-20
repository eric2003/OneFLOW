"""Join declared GoogleTest evidence to a configured build's CTest list.

This module performs CTest's non-executing list query only. Discovery is not
execution evidence and never changes a Contract's declared status.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys

from .contract_catalog import load_contract_catalog, resolve_contract_evidence
from .ctest_discovery import discover_configured_ctest, CTestDiscoveryResult, CTestTest
from .cmake_test_catalog import CMakeTestTarget, discover_cmake_test_targets


@dataclass(frozen=True)
class CTestEvidenceMatch:
    """Configured-build status for one declared GoogleTest evidence item."""

    contract_id: str
    reference: str
    source_path: str
    status: str
    targets: tuple[str, ...] = ()
    discovered_name: str = ""
    diagnostic: str = ""


def _target_provenance_matches(test: CTestTest, target: CMakeTestTarget) -> bool:
    """Require executable or CMake backtrace provenance for a target match."""
    if test.command:
        executable_name = Path(test.command[0].replace("\\", "/")).stem.casefold()
        if executable_name == target.target.casefold():
            return True
    if test.backtrace_file:
        backtrace = test.backtrace_file.replace("\\", "/").casefold()
        cmake_file = target.cmake_file.replace("\\", "/").casefold()
        if backtrace == cmake_file or backtrace.endswith("/" + cmake_file):
            return True
    return False


def match_configured_ctest_evidence(
    contract_id: str,
    reference: str,
    source_path: str,
    targets: tuple[CMakeTestTarget, ...],
    tests: tuple[CTestTest, ...],
) -> CTestEvidenceMatch:
    """Match one exact GoogleTest identity to one target-provenanced CTest item."""
    if not targets:
        return CTestEvidenceMatch(
            contract_id, reference, source_path, "target_unmapped",
            diagnostic="No supported CMake target maps this source file.",
        )
    declared = tuple(item for item in targets if item.discovery_declared)
    if not declared:
        return CTestEvidenceMatch(
            contract_id, reference, source_path, "discovery_not_declared",
            tuple(item.target for item in targets),
            diagnostic="Mapped CMake target has no gtest_discover_tests declaration.",
        )
    exact_named = tuple(item for item in tests if item.name == reference)
    if not exact_named:
        deferred = any(item.discovery_mode.casefold() == "pre_test" for item in declared)
        status = "discovery_deferred" if deferred else "not_discovered"
        diagnostic = (
            "This target uses PRE_TEST discovery; run CTest once, then query the list again."
            if deferred else
            "The exact GoogleTest identity is absent from this configured CTest list."
        )
        return CTestEvidenceMatch(
            contract_id, reference, source_path, status,
            tuple(item.target for item in declared), diagnostic=diagnostic,
        )
    matches = []
    for test in exact_named:
        for target in declared:
            if _target_provenance_matches(test, target):
                matches.append((test, target))
    if len(matches) == 1:
        test, target = matches[0]
        return CTestEvidenceMatch(
            contract_id, reference, source_path, "discovered",
            (target.target,), test.name,
            "Exact CTest name and CMake target provenance both match.",
        )
    if len(matches) > 1:
        return CTestEvidenceMatch(
            contract_id, reference, source_path, "ambiguous",
            tuple(sorted({target.target for _, target in matches})),
            diagnostic="More than one configured CTest entry/target matches this identity.",
        )
    return CTestEvidenceMatch(
        contract_id, reference, source_path, "target_provenance_mismatch",
        tuple(item.target for item in declared),
        diagnostic="CTest contains the exact name, but its executable/backtrace does not identify the mapped target.",
    )


def audit_configured_ctest_evidence(
    repository_root: Path, tests: tuple[CTestTest, ...]
) -> tuple[CTestEvidenceMatch, ...]:
    """Audit resolved GoogleTest contract evidence against configured CTest tests."""
    root = Path(repository_root).resolve()
    summaries, catalog_issues = load_contract_catalog(root)
    if catalog_issues:
        raise ValueError("Contract catalog errors: " + "; ".join(catalog_issues))
    results = []
    for summary in summaries:
        for resolution in resolve_contract_evidence(root, summary):
            if resolution.evidence.kind not in {"google_test", "linked_test"}:
                continue
            if resolution.status != "resolved":
                results.append(CTestEvidenceMatch(
                    summary.id, resolution.evidence.reference,
                    resolution.resolved_path, "stale_source",
                    diagnostic=resolution.diagnostic or "Contract source declaration is stale.",
                ))
                continue
            source = Path(resolution.resolved_path).resolve()
            try:
                relative = source.relative_to(root).as_posix()
            except ValueError:
                results.append(CTestEvidenceMatch(
                    summary.id, resolution.evidence.reference,
                    resolution.resolved_path, "stale_source",
                    diagnostic="Evidence source escapes the repository root.",
                ))
                continue
            targets = discover_cmake_test_targets(root, relative)
            results.append(match_configured_ctest_evidence(
                summary.id, resolution.evidence.reference, relative, targets, tests
            ))
    return tuple(results)


def format_ctest_evidence_report(
    discovery: CTestDiscoveryResult,
    matches: tuple[CTestEvidenceMatch, ...],
) -> str:
    """Render configured discovery separately from test execution outcomes."""
    lines = [
        "OneFLOW Configured CTest Evidence",
        f"Build directory: {discovery.build_directory}",
        f"Configuration: {discovery.configuration or '(generator default)'}",
        f"CTest version: {discovery.snapshot.version.text if discovery.snapshot else 'unknown'}",
        f"Configured CTest tests: {len(discovery.snapshot.tests) if discovery.snapshot else 0}",
        "No tests were executed by this command.",
    ]
    for item in matches:
        target_text = ", ".join(item.targets) if item.targets else "-"
        detail = f" -- {item.diagnostic}" if item.diagnostic else ""
        lines.append(
            f"[{item.status}] {item.contract_id} :: {item.reference} "
            f"(source={item.source_path or '-'}, target={target_text}){detail}"
        )
    counts: dict[str, int] = {}
    for item in matches:
        counts[item.status] = counts.get(item.status, 0) + 1
    lines.append("Summary: " + (", ".join(
        f"{name}={count}" for name, count in sorted(counts.items())
    ) if counts else "no GoogleTest evidence declarations"))
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Match GoogleTest contract evidence against a configured CTest list."
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--build-dir", type=Path, default=None)
    parser.add_argument("--configuration", default=None)
    parser.add_argument("--ctest", default="ctest")
    parser.add_argument("--require-all", action="store_true",
                        help="return 1 unless every declared GoogleTest is discovered")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    build = (args.build_dir or root / "build").expanduser().resolve()
    discovery = discover_configured_ctest(
        build, args.configuration, ctest_executable=args.ctest
    )
    if discovery.status != "discovered" or discovery.snapshot is None:
        print("OneFLOW Configured CTest Evidence")
        print(f"FAIL: {discovery.status}: {discovery.diagnostic}")
        return 2
    try:
        matches = audit_configured_ctest_evidence(root, discovery.snapshot.tests)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"OneFLOW Configured CTest Evidence\nFAIL: {exc}")
        return 2
    print(format_ctest_evidence_report(discovery, matches))
    incomplete = any(item.status != "discovered" for item in matches)
    return 1 if args.require_all and incomplete else 0


if __name__ == "__main__":
    sys.exit(main())
