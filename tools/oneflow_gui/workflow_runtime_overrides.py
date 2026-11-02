"""Apply runtime setting overrides to an isolated case parameter file."""

from __future__ import annotations

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
    overrides: Mapping[str, str],
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
        return f'{prefix}{name}{equals}"{_quote(overrides[name])}"{suffix}'

    content = _SETTING_PATTERN.sub(replace_setting, content)
    missing = [name for name in overrides if name not in updated]
    for name in missing:
        if not _IDENTIFIER_PATTERN.fullmatch(name):
            raise ValueError(f"Invalid runtime setting name: {name!r}")
        content += ("" if not content or content.endswith(("\n", "\r")) else newline)
        content += f'string {name} = "{_quote(overrides[name])}";{newline}'

    path.write_bytes(content.encode(encoding))
    return tuple(overrides)


def _quote(value: str) -> str:
    if not isinstance(value, str):
        raise TypeError("Runtime override values must be strings.")
    return value.replace("\\", "\\\\").replace('"', '\\"')
