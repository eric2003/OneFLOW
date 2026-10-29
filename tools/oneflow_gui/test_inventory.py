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


def _python_tests(repository_root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    tools_root = repository_root / "tools"
    if not tools_root.is_dir():
        return records
    for path in sorted(tools_root.rglob("test_*.py")):
        if "__pycache__" in path.parts or ".venv" in path.parts:
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (OSError, UnicodeDecodeError, SyntaxError):
            continue
        module_name = ".".join(path.relative_to(repository_root).with_suffix("").parts)
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for method in node.body:
                if isinstance(method, (ast.FunctionDef, ast.AsyncFunctionDef)) and method.name.startswith("test_"):
                    records.append({
                        "id": f"unittest:{module_name}.{node.name}.{method.name}",
                        "framework": "unittest",
                        "suite": node.name,
                        "name": method.name,
                        "kind": "test method",
                        "source": path.relative_to(repository_root).as_posix(),
                        "line": method.lineno,
                        "discovery": "source",
                        "scenarioId": None,
                        "scenarioRole": None,
                        "input": None,
                        "expected": None,
                        "result": None,
                        "execution": f"python -m unittest {module_name}.{node.name}.{method.name}",
                    })
    return records


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
