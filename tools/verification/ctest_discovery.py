"""Configured-build CTest discovery; this module never executes tests."""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import subprocess
from typing import Any


class CTestDiscoveryError(ValueError):
    """Invalid CTest JSON or unsupported schema."""


class CTestDiscoveryJSONError(CTestDiscoveryError):
    """CTest did not emit valid JSON."""


@dataclass(frozen=True)
class CTestVersion:
    major: int
    minor: int
    patch: int

    @property
    def text(self) -> str:
        return f"{self.major}.{self.minor}.{self.patch}"


@dataclass(frozen=True)
class CTestTest:
    name: str
    configuration: str | None
    command: tuple[str, ...] | None
    properties: tuple[dict[str, Any], ...]
    backtrace_file: str | None = None
    backtrace_line: int | None = None


@dataclass(frozen=True)
class CTestDiscoverySnapshot:
    version: CTestVersion
    tests: tuple[CTestTest, ...]


@dataclass(frozen=True)
class CTestDiscoveryResult:
    status: str
    build_directory: str
    configuration: str | None
    snapshot: CTestDiscoverySnapshot | None = None
    command: tuple[str, ...] = ()
    return_code: int | None = None
    stdout: str = ""
    stderr: str = ""
    diagnostic: str = ""


def _integer(value: Any, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise CTestDiscoveryError(f"{label} must be a non-negative integer")
    return value


def _backtrace(test: dict[str, Any], graph: dict[str, Any] | None):
    if "backtrace" not in test or graph is None:
        return None, None
    nodes, files = graph.get("nodes"), graph.get("files")
    if not isinstance(nodes, list) or not isinstance(files, list):
        raise CTestDiscoveryError("backtraceGraph must contain nodes and files arrays")
    index = _integer(test["backtrace"], "test.backtrace")
    source_file = line = None
    visited: set[int] = set()
    while True:
        if index >= len(nodes):
            raise CTestDiscoveryError("backtraceGraph references a missing node")
        if index in visited:
            raise CTestDiscoveryError("backtraceGraph contains a parent cycle")
        visited.add(index)
        node = nodes[index]
        if not isinstance(node, dict):
            raise CTestDiscoveryError("backtraceGraph.nodes entries must be objects")
        file_index = node.get("file")
        if source_file is None and file_index is not None:
            file_index = _integer(file_index, "backtrace node file")
            if file_index >= len(files) or not isinstance(files[file_index], str):
                raise CTestDiscoveryError("backtrace node references an invalid file")
            source_file = files[file_index]
        node_line = node.get("line")
        if line is None and node_line is not None:
            line = _integer(node_line, "backtrace node line")
            if line == 0:
                raise CTestDiscoveryError("backtrace node line must be positive")
        if source_file is not None and line is not None:
            break
        parent = node.get("parent")
        if parent is None:
            break
        index = _integer(parent, "backtrace node parent")
    return source_file, line


def parse_ctest_discovery_json(source: str) -> CTestDiscoverySnapshot:
    """Parse CTest's ctestInfo JSON without inferring execution outcomes."""
    try:
        payload = json.loads(source)
    except json.JSONDecodeError as exc:
        raise CTestDiscoveryJSONError(f"invalid CTest JSON: {exc}") from exc
    if not isinstance(payload, dict) or payload.get("kind") != "ctestInfo":
        raise CTestDiscoveryError("CTest JSON root kind must be 'ctestInfo'")
    version_data = payload.get("version")
    if not isinstance(version_data, dict):
        raise CTestDiscoveryError("CTest JSON version must be an object")
    version = CTestVersion(
        _integer(version_data.get("major"), "version.major"),
        _integer(version_data.get("minor"), "version.minor"),
        _integer(version_data.get("patch"), "version.patch"),
    )
    tests_data = payload.get("tests")
    if not isinstance(tests_data, list):
        raise CTestDiscoveryError("CTest JSON tests must be an array")
    graph = payload.get("backtraceGraph")
    if graph is not None and not isinstance(graph, dict):
        raise CTestDiscoveryError("backtraceGraph must be an object")

    tests = []
    for position, item in enumerate(tests_data):
        if not isinstance(item, dict):
            raise CTestDiscoveryError(f"tests[{position}] must be an object")
        name = item.get("name")
        if not isinstance(name, str) or not name.strip():
            raise CTestDiscoveryError(f"tests[{position}].name must be a non-empty string")
        config = item.get("config")
        if config is not None and not isinstance(config, str):
            raise CTestDiscoveryError(f"tests[{position}].config must be a string")
        raw_command = item.get("command")
        if raw_command is None:
            command = None
        elif isinstance(raw_command, list) and all(isinstance(x, str) for x in raw_command):
            command = tuple(raw_command)
        else:
            raise CTestDiscoveryError(f"tests[{position}].command must be an array of strings")
        properties = item.get("properties", [])
        if not isinstance(properties, list) or not all(
            isinstance(x, dict) and isinstance(x.get("name"), str) for x in properties
        ):
            raise CTestDiscoveryError(f"tests[{position}].properties must contain named objects")
        source_file, source_line = _backtrace(item, graph)
        tests.append(CTestTest(
            name=name, configuration=config, command=command,
            properties=tuple(properties), backtrace_file=source_file,
            backtrace_line=source_line,
        ))
    return CTestDiscoverySnapshot(version=version, tests=tuple(tests))


def discover_configured_ctest(
    build_directory: Path,
    configuration: str | None = None,
    *,
    ctest_executable: str = "ctest",
    timeout_seconds: float = 30.0,
) -> CTestDiscoveryResult:
    """Query a configured build tree using CTest's non-executing list mode."""
    build_dir = Path(build_directory).expanduser().resolve()
    if not build_dir.is_dir():
        return CTestDiscoveryResult(
            "build_directory_missing", str(build_dir), configuration,
            diagnostic="The selected build directory does not exist.",
        )
    executable = shutil.which(ctest_executable)
    if executable is None:
        return CTestDiscoveryResult(
            "ctest_unavailable", str(build_dir), configuration,
            diagnostic=f"CTest executable not found: {ctest_executable}",
        )
    command = [executable]
    if configuration:
        command.extend(["-C", configuration])
    command.extend(["-N", "--show-only=json-v1"])
    try:
        completed = subprocess.run(
            command, cwd=build_dir, capture_output=True, text=True,
            encoding="utf-8", errors="replace", timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        return CTestDiscoveryResult(
            "timeout", str(build_dir), configuration, command=tuple(command),
            stdout=exc.stdout or "", stderr=exc.stderr or "",
            diagnostic=f"CTest discovery exceeded {timeout_seconds:g} seconds.",
        )
    except OSError as exc:
        return CTestDiscoveryResult(
            "ctest_unavailable", str(build_dir), configuration,
            command=tuple(command), diagnostic=str(exc),
        )
    common = dict(
        build_directory=str(build_dir), configuration=configuration,
        command=tuple(command), return_code=completed.returncode,
        stdout=completed.stdout, stderr=completed.stderr,
    )
    if completed.returncode:
        return CTestDiscoveryResult(
            "command_failed", diagnostic="CTest returned a non-zero exit code.", **common
        )
    try:
        snapshot = parse_ctest_discovery_json(completed.stdout)
    except CTestDiscoveryJSONError as exc:
        return CTestDiscoveryResult("invalid_json", diagnostic=str(exc), **common)
    except CTestDiscoveryError as exc:
        return CTestDiscoveryResult("invalid_schema", diagnostic=str(exc), **common)
    return CTestDiscoveryResult(
        "discovered", snapshot=snapshot,
        diagnostic=f"Discovered {len(snapshot.tests)} CTest test(s); none were executed.",
        **common,
    )
