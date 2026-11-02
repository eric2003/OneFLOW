"""Discover repository tests and normalize metadata for Workbench views.

Source code remains the authority for executable tests. Scenario JSON adds
input/expected examples only when it explicitly links a test.
"""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path
from typing import Any

CPP_TEST_PATTERN = re.compile(
    r"\b(TEST|TEST_F|TEST_P)\s*\(\s*([A-Za-z_][A-Za-z_0-9]*)\s*,"
    r"\s*([A-Za-z_][A-Za-z_0-9]*)\s*\)"
)


def _scenario_links(repository_root: Path) -> dict[tuple[str, str, str], dict[str, Any]]:
    links: dict[tuple[str, str, str], dict[str, Any]] = {}
    scenario_root = repository_root / "tools" / "geometry_workbench" / "scenarios"
    if not scenario_root.is_dir():
        return links
    for path in sorted(scenario_root.glob("*.json")):
        try:
            scenario = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for role, key in (("contract", "test"), ("integration", "integrationTest")):
            test = scenario.get(key) or {}
            identity = (str(test.get("framework", "GoogleTest")),
                        str(test.get("suite", "")), str(test.get("name", "")))
            if identity[1] and identity[2]:
                links[identity] = {
                    "scenarioId": str(scenario.get("id", path.stem)),
                    "scenarioRole": role,
                    "input": scenario.get("input") if role == "contract" else None,
                    "expected": scenario.get("expected") if role == "contract" else None,
                }
    return links


def _cpp_tests(repository_root: Path, links: dict) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    test_root = repository_root / "tests"
    if not test_root.is_dir():
        return records
    for path in sorted(test_root.rglob("*.cpp")):
        try:
            source = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for match in CPP_TEST_PATTERN.finditer(source):
            macro, suite, name = match.groups()
            identity = ("GoogleTest", suite, name)
            link = links.get(identity, {})
            records.append({
                "id": f"gtest:{suite}.{name}",
                "framework": "GoogleTest",
                "suite": suite,
                "name": name,
                "kind": macro,
                "source": path.relative_to(repository_root).as_posix(),
                "line": source.count("\n", 0, match.start()) + 1,
                "discovery": "source",
                "scenarioId": link.get("scenarioId"),
                "scenarioRole": link.get("scenarioRole"),
                "input": link.get("input"),
                "expected": link.get("expected"),
                "result": None,
                "execution": "CTest target; per-test filtering depends on target registration",
            })
    return records


def _python_test_roots(repository_root: Path) -> tuple[Path, ...]:
    return tuple(
        repository_root / relative
        for relative in ("tools", "tests", Path("gui/python/tkinter"))
        if (repository_root / relative).is_dir()
    )


def _has_python_test_declarations(tree: ast.AST) -> bool:
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            if any(
                isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef))
                and method.name.startswith("test_")
                for method in node.body
            ):
                return True
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("test_"):
                return True
    return False


def _python_tests(repository_root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for source_root in _python_test_roots(repository_root):
        for path in sorted(source_root.rglob("test_*.py")):
            if "__pycache__" in path.parts or ".venv" in path.parts:
                continue
            try:
                source = path.read_text(encoding="utf-8")
                tree = ast.parse(source, filename=str(path))
            except (OSError, UnicodeDecodeError, SyntaxError):
                continue
            module_name = ".".join(path.relative_to(repository_root).with_suffix("").parts)
            relative_source = path.relative_to(repository_root).as_posix()
            for node in tree.body:
                if isinstance(node, ast.ClassDef):
                    for method in node.body:
                        if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef)) and method.name.startswith("test_"):
                            records.append({
                                "id": f"unittest:{module_name}.{node.name}.{method.name}",
                                "framework": "unittest",
                                "suite": node.name,
                                "name": method.name,
                                "kind": "test method",
                                "source": relative_source,
                                "line": method.lineno,
                                "discovery": "source",
                                "scenarioId": None,
                                "scenarioRole": None,
                                "input": None,
                                "expected": None,
                                "result": None,
                                "execution": f"python -m unittest {module_name}.{node.name}.{method.name}",
                            })
                elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
                    records.append({
                        "id": f"unittest:{module_name}.{node.name}",
                        "framework": "unittest",
                        "suite": module_name,
                        "name": node.name,
                        "kind": "test function",
                        "source": relative_source,
                        "line": node.lineno,
                        "discovery": "source",
                        "scenarioId": None,
                        "scenarioRole": None,
                        "input": None,
                        "expected": None,
                        "result": None,
                        "execution": f"python -m unittest {module_name}.{node.name}",
                    })
    return records


def audit_python_test_registry(
    repository_root: Path,
    commands: tuple | None = None,
    exclusions: dict[str, str] | None = None,
) -> dict[str, tuple[str, ...]]:
    """Compare discovered unittest source files with the explicit execution registry."""
    from .test_runner import _PYTHON_TEST_EXCLUSIONS, _PYTHON_UNIT_TESTS

    root = repository_root.resolve()
    command_entries = _PYTHON_UNIT_TESTS if commands is None else commands
    exclusion_reasons = _PYTHON_TEST_EXCLUSIONS if exclusions is None else exclusions

    discovered: set[str] = set()
    for source_root in _python_test_roots(root):
        for path in source_root.rglob("test_*.py"):
            try:
                tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            except (OSError, UnicodeDecodeError, SyntaxError):
                continue
            if _has_python_test_declarations(tree):
                discovered.add(path.relative_to(root).as_posix())

    registered: set[str] = set()
    for _label, arguments in command_entries:
        if len(arguments) >= 3 and arguments[0:2] == ("-m", "unittest"):
            if arguments[2] == "discover":
                start_dir = Path(".")
                pattern = "test*.py"
                for index, argument in enumerate(arguments[:-1]):
                    if argument == "-s":
                        start_dir = Path(arguments[index + 1])
                    elif argument == "-p":
                        pattern = arguments[index + 1]
                search_root = (root / start_dir).resolve()
                if search_root.is_dir():
                    registered.update(
                        path.relative_to(root).as_posix()
                        for path in search_root.rglob(pattern)
                        if path.is_file()
                    )
            else:
                module_path = root.joinpath(*arguments[2].split(".")).with_suffix(".py")
                if module_path.is_file():
                    registered.add(module_path.relative_to(root).as_posix())
        elif arguments:
            script_path = (root / arguments[0]).resolve()
            if script_path.is_file():
                try:
                    registered.add(script_path.relative_to(root).as_posix())
                except ValueError:
                    pass

    excluded = discovered.intersection(exclusion_reasons)
    unregistered = discovered.difference(registered).difference(excluded)
    stale_exclusions = set(exclusion_reasons).difference(discovered)
    return {
        "discovered": tuple(sorted(discovered)),
        "registered": tuple(sorted(discovered.intersection(registered))),
        "excluded": tuple(sorted(excluded)),
        "unregistered": tuple(sorted(unregistered)),
        "stale_exclusions": tuple(sorted(stale_exclusions)),
    }

def discover_tests(repository_root: Path) -> list[dict[str, Any]]:
    """Return a stable catalog of C++ and Python tests.

    Input, expected and result are separate fields. Discovery never invents
    values for tests that have no data-driven contract.
    """
    root = repository_root.resolve()
    links = _scenario_links(root)
    records = _cpp_tests(root, links) + _python_tests(root)
    return sorted(records, key=lambda item: (
        item["framework"], item["source"], item["line"], item["name"]
    ))
