"""Static mapping from GoogleTest source files to CMake test targets.

This module reads repository CMake declarations only. It does not configure,
compile, inspect a build directory, run CTest, or infer test outcomes.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shlex


@dataclass(frozen=True)
class CMakeTestTarget:
    """A test target whose declared source list contains a GoogleTest file."""

    target: str
    cmake_file: str
    source_file: str
    discovery_registered: bool
    discovery_mode: str = ""


def _cmake_commands(source: str):
    """Yield command name and raw argument text, ignoring comments."""
    index = 0
    length = len(source)
    command_start = re.compile(r"[A-Za-z_][A-Za-z_0-9]*")

    while index < length:
        if source[index] == "#":
            newline = source.find("\n", index)
            index = length if newline < 0 else newline + 1
            continue
        match = command_start.match(source, index)
        if match is None:
            index += 1
            continue
        name = match.group()
        cursor = match.end()
        while cursor < length and source[cursor].isspace():
            cursor += 1
        if cursor >= length or source[cursor] != "(":
            index = match.end()
            continue

        start = cursor + 1
        cursor += 1
        depth = 1
        quote = False
        escaped = False
        while cursor < length and depth:
            current = source[cursor]
            if escaped:
                escaped = False
            elif current == "\\" and quote:
                escaped = True
            elif current == '"':
                quote = not quote
            elif not quote and current == "#":
                newline = source.find("\n", cursor)
                cursor = length if newline < 0 else newline
                continue
            elif not quote and current == "(":
                depth += 1
            elif not quote and current == ")":
                depth -= 1
            cursor += 1
        if depth == 0:
            yield name.lower(), source[start:cursor - 1]
            index = cursor
        else:
            break


def _arguments(raw: str) -> list[str]:
    try:
        return shlex.split(raw, comments=True, posix=True)
    except ValueError:
        return []


def _source_path(token: str, cmake_directory: Path, repository_root: Path) -> Path | None:
    normalized = token.replace("\\", "/")
    for prefix in ("${CMAKE_SOURCE_DIR}/", "${PROJECT_SOURCE_DIR}/"):
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix):]
            return (repository_root / normalized).resolve()
    if "${" in normalized or "$<" in normalized:
        return None
    candidate = Path(normalized)
    if candidate.suffix.lower() not in {".cpp", ".cc", ".cxx"}:
        return None
    if not candidate.is_absolute():
        candidate = cmake_directory / candidate
    try:
        return candidate.resolve()
    except OSError:
        return None


def discover_cmake_test_targets(
    repository_root: Path,
    source_file: str,
) -> tuple[CMakeTestTarget, ...]:
    """Find direct add_executable source mappings and CTest discovery declarations.

    Sources supplied through variables, target_sources(), functions, or
    generated lists are intentionally not guessed by this static reader.
    """
    root = repository_root.resolve()
    requested = Path(source_file)
    if requested.is_absolute():
        try:
            requested = requested.resolve().relative_to(root)
        except ValueError:
            return ()
    requested_path = (root / requested).resolve()

    targets: dict[str, tuple[Path, set[Path]]] = {}
    discovery: dict[str, str] = {}
    tests_root = root / "tests"
    if not tests_root.is_dir():
        return ()

    for cmake_path in sorted(tests_root.rglob("CMakeLists.txt")):
        try:
            content = cmake_path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for command, raw in _cmake_commands(content):
            arguments = _arguments(raw)
            if command == "add_executable" and len(arguments) >= 2:
                target = arguments[0]
                sources = {
                    resolved
                    for token in arguments[1:]
                    if (resolved := _source_path(token, cmake_path.parent, root)) is not None
                }
                if sources:
                    previous = targets.get(target)
                    if previous is None:
                        targets[target] = (cmake_path, sources)
                    else:
                        previous[1].update(sources)
            elif command == "gtest_discover_tests" and arguments:
                mode = ""
                for index, argument in enumerate(arguments[:-1]):
                    if argument == "DISCOVERY_MODE":
                        mode = arguments[index + 1]
                        break
                discovery[arguments[0]] = mode

    matches = []
    for target, (cmake_path, sources) in sorted(targets.items()):
        if requested_path not in sources:
            continue
        try:
            cmake_relative = cmake_path.relative_to(root).as_posix()
            source_relative = requested_path.relative_to(root).as_posix()
        except ValueError:
            continue
        mode = discovery.get(target, "")
        matches.append(
            CMakeTestTarget(
                target=target,
                cmake_file=cmake_relative,
                source_file=source_relative,
                discovery_registered=target in discovery,
                discovery_mode=mode,
            )
        )
    return tuple(matches)
