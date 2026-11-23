"""Repository-wide consistency checks to run before committing changes.

The preflight is read-only. It reuses the Contract Catalog and Python test
inventory instead of maintaining a second list of contracts or test modules.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable

from .contract_catalog import (
    _has_python_test_declaration,
    _mask_cpp_non_code,
    _read_python_source,
    load_contract_catalog,
    resolve_contract_evidence,
)
from tools.oneflow_gui.test_inventory import audit_python_test_registry
from .cmake_test_catalog import discover_cmake_test_targets


def _changed_files(root: Path) -> tuple[str, ...]:
    """Return tracked, staged, and untracked paths relative to the repository."""
    commands = (
        ("git", "diff", "--name-only", "HEAD"),
        ("git", "diff", "--cached", "--name-only"),
        ("git", "ls-files", "--others", "--exclude-standard"),
    )
    paths: set[str] = set()
    for command in commands:
        result = subprocess.run(
            command, cwd=root, text=True, capture_output=True, check=False
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or "git command failed"
            raise RuntimeError(f"cannot determine changed files: {detail}")
        paths.update(line.strip().replace("\\", "/")
                     for line in result.stdout.splitlines() if line.strip())
    return tuple(sorted(paths))



def find_python_test_candidates(
    repository_root: Path, reference: str
) -> tuple[str, ...]:
    """Find repository Python files declaring an exact unittest identity.

    This is diagnostic only: candidates are never written back to a Contract.
    """
    root = Path(repository_root).resolve()
    ignored = {".git", ".venv", "venv", "build", "__pycache__", ".pytest_cache"}
    matches: list[str] = []
    for path in root.rglob("*.py"):
        try:
            relative = path.relative_to(root)
            if any(part in ignored for part in relative.parts):
                continue
            resolved = path.resolve()
            resolved.relative_to(root)
            if not resolved.is_file():
                continue
            source = _read_python_source(resolved)
            if _has_python_test_declaration(source, reference):
                matches.append(relative.as_posix())
        except (OSError, UnicodeDecodeError, SyntaxError, LookupError, ValueError):
            continue
    return tuple(sorted(matches))


def _stale_python_evidence_diagnostic(
    root: Path, reference: str, original_diagnostic: str
) -> str:
    """Add exact declaration-based migration hints to a stale evidence error."""
    candidates = find_python_test_candidates(root, reference)
    if len(candidates) == 1:
        return (
            f"{original_diagnostic}; matching declaration found at "
            f"{candidates[0]} (update evidence scope after review)"
        )
    if len(candidates) > 1:
        return (
            f"{original_diagnostic}; matching declarations found at: "
            f"{', '.join(candidates)} (disambiguate evidence scope)"
        )
    return original_diagnostic



def find_cpp_test_candidates(
    repository_root: Path, reference: str
) -> tuple[str, ...]:
    """Find C++ files containing an exact GoogleTest declaration."""
    suite, separator, test_name = reference.partition(".")
    if not separator or not suite or not test_name:
        return ()
    declaration = re.compile(
        rf"\b(?:TEST(?:_F|_P)?|TYPED_TEST(?:_P)?)\s*"
        rf"\(\s*{re.escape(suite)}\s*,\s*{re.escape(test_name)}\s*\)"
    )
    root = Path(repository_root).resolve()
    ignored = {".git", ".venv", "venv", "build", "out", "cmake-build-debug",
               "cmake-build-release", "__pycache__"}
    matches: list[str] = []
    for suffix in ("*.cpp", "*.cc", "*.cxx"):
        for path in root.rglob(suffix):
            try:
                relative = path.relative_to(root)
                if any(part.lower() in ignored for part in relative.parts):
                    continue
                resolved = path.resolve()
                resolved.relative_to(root)
                if not resolved.is_file():
                    continue
                source = _mask_cpp_non_code(
                    resolved.read_text(encoding="utf-8")
                )
                if declaration.search(source):
                    matches.append(relative.as_posix())
            except (OSError, UnicodeDecodeError, ValueError):
                continue
    return tuple(sorted(set(matches)))


def _stale_cpp_evidence_diagnostic(
    root: Path, reference: str, original_diagnostic: str
) -> str:
    """Add exact GoogleTest declaration migration hints to a stale evidence error."""
    candidates = find_cpp_test_candidates(root, reference)
    if len(candidates) == 1:
        return (
            f"{original_diagnostic}; matching GoogleTest declaration found at "
            f"{candidates[0]} (update evidence scope after review)"
        )
    if len(candidates) > 1:
        return (
            f"{original_diagnostic}; matching GoogleTest declarations found at: "
            f"{', '.join(candidates)} (disambiguate evidence scope)"
        )
    return original_diagnostic


def _cpp_test_target_diagnostic(root: Path, source_path: str) -> str:
    """Check discovery for a statically mapped GoogleTest source."""
    records = discover_cmake_test_targets(root, source_path)
    missing = tuple(record for record in records if not record.discovery_declared)
    if not missing:
        return ""
    details = ", ".join(f"{record.target} ({record.cmake_file})" for record in missing)
    return f"CMake test target mapping found without gtest_discover_tests: {details}"


def collect_preflight_issues(repository_root: Path) -> tuple[str, ...]:
    """Validate repository contracts, evidence declarations, and test coverage."""
    root = Path(repository_root).resolve()
    issues: list[str] = []
    summaries, catalog_issues = load_contract_catalog(root)
    issues.extend(f"catalog: {item}" for item in catalog_issues)

    registry = audit_python_test_registry(root)
    for path in registry["unregistered"]:
        issues.append(
            f"python test is not covered by the execution registry: {path}"
        )
    for path in registry["stale_exclusions"]:
        issues.append(f"test exclusion no longer matches a discovered test: {path}")

    for summary in summaries:
        seen: set[tuple[str, str]] = set()
        for evidence in summary.evidence:
            identity = (evidence.kind, evidence.reference)
            if identity in seen:
                issues.append(
                    f"{summary.id}: duplicate evidence declaration "
                    f"{evidence.kind}: {evidence.reference}"
                )
            seen.add(identity)

        for resolution in resolve_contract_evidence(root, summary):
            if resolution.status == "stale":
                diagnostic = resolution.diagnostic
                if resolution.evidence.kind == "python_test":
                    diagnostic = _stale_python_evidence_diagnostic(
                        root, resolution.evidence.reference, diagnostic
                    )
                elif resolution.evidence.kind in ("linked_test", "google_test"):
                    diagnostic = _stale_cpp_evidence_diagnostic(
                        root, resolution.evidence.reference, diagnostic
                    )
                issues.append(
                    f"{summary.id}: stale {resolution.evidence.kind} evidence "
                    f"{resolution.evidence.reference}: {diagnostic}"
                )
                continue
            if resolution.evidence.kind == "python_test" and resolution.status == "resolved":
                source = Path(resolution.resolved_path).resolve()
                try:
                    relative = source.relative_to(root).as_posix()
                except ValueError:
                    issues.append(
                        f"{summary.id}: Python evidence source escapes repository: "
                        f"{resolution.resolved_path}"
                    )
                    continue
                if relative not in registry["registered"]:
                    issues.append(
                        f"{summary.id}: Python evidence is not covered by the "
                        f"Python execution registry: {relative} "
                        f"({resolution.evidence.reference})"
                    )
            elif resolution.status == "resolved" and resolution.evidence.kind in ("linked_test", "google_test"):
                source = Path(resolution.resolved_path).resolve()
                try:
                    relative = source.relative_to(root).as_posix()
                except ValueError:
                    issues.append(f"{summary.id}: GoogleTest evidence source escapes repository: {resolution.resolved_path}")
                    continue
                diagnostic = _cpp_test_target_diagnostic(root, relative)
                if diagnostic:
                    issues.append(f"{summary.id}: {resolution.evidence.kind} evidence {resolution.evidence.reference}: {diagnostic}")

    return tuple(dict.fromkeys(issues))


def format_report(
    issues: Iterable[str], changed_files: Iterable[str] = ()
) -> str:
    """Render a concise report suitable for terminals and CI logs."""
    issue_list = tuple(issues)
    changed = tuple(changed_files)
    lines = ["OneFLOW Change Preflight"]
    if changed:
        lines.append(f"Changed files: {len(changed)}")
        lines.extend(f"  - {path}" for path in changed)
    if issue_list:
        lines.append(f"FAIL: {len(issue_list)} issue(s)")
        lines.extend(f"  - {issue}" for issue in issue_list)
    else:
        lines.append("PASS: contract catalog, evidence references, and test registries are consistent.")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check OneFLOW contract and test-registration consistency."
    )
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[2],
        help="repository root (defaults to the current OneFLOW checkout)",
    )
    parser.add_argument(
        "--changed", action="store_true",
        help="also list changed paths; consistency checks still cover the full catalog",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()
    try:
        changed = _changed_files(root) if args.changed else ()
        issues = collect_preflight_issues(root)
    except (OSError, RuntimeError) as exc:
        print(f"OneFLOW Change Preflight\nFAIL: {exc}")
        return 2
    print(format_report(issues, changed))
    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
