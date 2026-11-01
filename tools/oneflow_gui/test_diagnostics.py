"""Structured diagnostics parsed from Python unittest and GoogleTest output.

This module does not execute tests or infer contract links. Callers supply the
actual command context and any explicitly resolved contract identifiers.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True, slots=True)
class TestDiagnostic:
    """One failed test with its execution context and optional contract links."""

    task_label: str
    framework: str
    test_id: str | None
    test_name: str | None
    source_file: str | None
    source_line: int | None
    exception_type: str | None
    message: str | None
    configuration: str | None
    working_directory: str | None
    started_at: str | None
    finished_at: str | None
    exit_code: int | None
    contract_ids: tuple[str, ...] = ()
    verification_impact: str = "unlinked"

    def to_dict(self) -> dict:
        """Return a JSON-compatible representation with stable field names."""
        payload = asdict(self)
        payload["contract_ids"] = list(self.contract_ids)
        return payload


def _source_location(lines: list[str]) -> tuple[str | None, int | None]:
    for line in lines:
        match = re.search(r'File "([^"]+)", line (\d+)', line)
        if match:
            return match.group(1), int(match.group(2))
        match = re.search(r"((?:[A-Za-z]:)?[^\s:]+\.(?:cpp|cc|cxx|h|hpp)):(\d+)(?::\d+)?", line)
        if match:
            return match.group(1), int(match.group(2))
    return None, None


def _exception(lines: list[str]) -> tuple[str | None, str | None]:
    for line in reversed(lines):
        match = re.match(r"^\s*([A-Za-z_][A-Za-z_0-9]*(?:Error|Exception)):\s*(.*)$", line)
        if match:
            return match.group(1), match.group(2).strip() or None
    for line in reversed(lines):
        match = re.match(r"^\s*([A-Za-z_][A-Za-z_0-9]*(?:Error|Exception))\s*$", line)
        if match:
            return match.group(1), None
    return None, None


def parse_test_diagnostics(
    output_lines: Iterable[str],
    *,
    task_label: str,
    framework: str,
    configuration: str | None = None,
    working_directory: str | Path | None = None,
    started_at: str | None = None,
    finished_at: str | None = None,
    exit_code: int | None = None,
    contract_ids: Iterable[str] = (),
) -> tuple[TestDiagnostic, ...]:
    """Parse failed unittest or GoogleTest cases without guessing missing data."""
    lines = list(output_lines)
    contracts = tuple(dict.fromkeys(str(item) for item in contract_ids if str(item)))
    impact = "requires_review" if contracts else "unlinked"
    records: list[TestDiagnostic] = []

    if framework.lower() in {"unittest", "python unittest"}:
        starts = [
            index for index, line in enumerate(lines)
            if re.match(r"^\s*(?:FAIL|ERROR):\s+", line)
        ]
        for position, start in enumerate(starts):
            end = starts[position + 1] if position + 1 < len(starts) else len(lines)
            block = lines[start:end]
            header = re.match(r"^\s*(FAIL|ERROR):\s+(.+?)\s*$", block[0])
            if header is None:
                continue
            name = header.group(2)
            source, line_number = _source_location(block)
            error_type, message = _exception(block)
            records.append(TestDiagnostic(
                task_label, "unittest", name, name, source, line_number,
                error_type, message, configuration,
                str(Path(working_directory).resolve()) if working_directory else None,
                started_at, finished_at, exit_code, contracts, impact,
            ))

    elif framework.lower() in {"googletest", "google test", "gtest"}:
        for index, line in enumerate(lines):
            match = re.match(r"^\s*\[\s*FAILED\s*\]\s+(.+?)(?:\s+\([^)]*\))?\s*$", line)
            if not match:
                continue
            name = match.group(1)
            block_start = max(0, index - 12)
            block = lines[block_start:index + 1]
            source, line_number = _source_location(block)
            error_type, message = _exception(block)
            if error_type is None:
                reason_lines = [
                    item.strip() for item in block
                    if item.strip()
                    and not item.strip().startswith("[ RUN")
                    and not item.strip().startswith("[ FAILED")
                    and not re.search(r"((?:[A-Za-z]:)?[^\s:]+\.(?:cpp|cc|cxx|h|hpp)):(\d+)", item)
                ]
                if reason_lines:
                    error_type = "GoogleTestFailure"
                    message = " | ".join(reason_lines[-4:])
            records.append(TestDiagnostic(
                task_label, "GoogleTest", name, name, source, line_number,
                error_type, message, configuration,
                str(Path(working_directory).resolve()) if working_directory else None,
                started_at, finished_at, exit_code, contracts, impact,
            ))

    return tuple(records)
