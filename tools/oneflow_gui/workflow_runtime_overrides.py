"""Apply runtime setting overrides to an isolated case parameter file."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Mapping


_SETTING_PATTERN = re.compile(
    r"^(\s*(?:string|int|Real|real|bool)\s+)"
    r"([A-Za-z_][A-Za-z0-9_]*)"
    r"(\s*=\s*)(.*?)"
    r"(;\s*(?://.*)?)$",
    re.MULTILINE,
)
_IDENTIFIER_PATTERN = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def apply_runtime_overrides(
    parameter_file: Path | str,
    overrides: Mapping[str, str | int | float | bool],
) -> tuple[str, ...]:
    """Update settings in a run-specific parameter file, never the source case.

    Existing declarations keep their type, indentation, and trailing comments.
    Missing settings are appended as string declarations.
    """
    path = Path(parameter_file)
    if not path.is_file():
        raise FileNotFoundError(f"Parameter file not found: {path}")
    if not overrides:
        return ()

    raw_content = path.read_bytes()
    if raw_content.startswith(b"\xef\xbb\xbf"):
        encoding = "utf-8-sig"
    else:
        try:
            raw_content.decode("utf-8")
            encoding = "utf-8"
        except UnicodeDecodeError:
            encoding = "gbk"
    content = raw_content.decode(encoding)

    newline = "\r\n" if "\r\n" in content else "\n"
    updated: set[str] = set()

    def replace_setting(match: re.Match[str]) -> str:
        prefix, name, equals, value, suffix = match.groups()
        if name not in overrides:
            return match.group(0)
        updated.add(name)
        declaration = prefix.strip().split()[0]
        rendered = _format_override(overrides[name], declaration)
        return f"{prefix}{name}{equals}{rendered}{suffix}"

    content = _SETTING_PATTERN.sub(replace_setting, content)
    missing = [name for name in overrides if name not in updated]
    for name in missing:
        if not _IDENTIFIER_PATTERN.fullmatch(name):
            raise ValueError(f"Invalid runtime setting name: {name!r}")
        content += ("" if not content or content.endswith(("\n", "\r")) else newline)
        value = overrides[name]
        if not isinstance(value, str):
            raise ValueError(f"Missing setting {name!r} can only be appended as a string.")
        content += f'string {name} = "{_quote(value)}";{newline}'

    path.write_bytes(content.encode(encoding))
    return tuple(overrides)


def _quote(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("Runtime override values must be strings.")
    return value.replace("\\", "\\\\").replace('"', '\\"')


def _format_override(value: str | int | float | bool, declaration: str) -> str:
    """Render a typed runtime value without changing the parameter declaration."""
    if declaration == "string":
        if not isinstance(value, str):
            raise ValueError("String parameter overrides require string values.")
        return f'"{_quote(value)}"'
    if declaration == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError("Integer parameter overrides require integer values.")
        return str(value)
    if declaration in {"Real", "real"}:
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise ValueError("Real parameter overrides require finite numeric values.")
        return format(float(value), ".15g")
    if declaration == "bool":
        if not isinstance(value, bool):
            raise ValueError("Boolean parameter overrides require boolean values.")
        return "true" if value else "false"
    raise ValueError(f"Unsupported parameter declaration type: {declaration!r}")
