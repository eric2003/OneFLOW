"""Lossless parser and writer for legacy OneFLOW parameter files."""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

_DECL = re.compile(
    r"^[ \t]*(?P<type>int|real|string)[ \t]+(?P<name>[A-Za-z_][A-Za-z0-9_]*)"
    r"[ \t]*=[ \t]*(?P<value>\"(?:\\.|[^\"\\])*\"|'(?:\\.|[^'\\])*'|"
    r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)"
    r"(?P<tail>[ \t]*(?:;)?[ \t]*(?://[^\r\n]*)?)[ \t]*\r?$",
    re.MULTILINE | re.IGNORECASE,
)


@dataclass(frozen=True, slots=True)
class ParameterDeclaration:
    name: str
    declared_type: str
    raw_value: str
    value: Any
    value_start: int
    value_end: int


@dataclass(frozen=True, slots=True)
class ParameterChoice:
    label: str
    value: Any


@dataclass(frozen=True, slots=True)
class ParameterSpec:
    key: str
    label: str
    description: str
    value_kind: str
    group: str
    choices: tuple[ParameterChoice, ...] = ()
    minimum: float | None = None
    maximum: float | None = None
    decimals: int = 6


def _decode(data: bytes) -> tuple[str, str]:
    if data.startswith(b"\xef\xbb\xbf"):
        return data.decode("utf-8-sig"), "utf-8-sig"
    try:
        return data.decode("utf-8"), "utf-8"
    except UnicodeDecodeError:
        try:
            return data.decode("gbk"), "gbk"
        except UnicodeDecodeError as error:
            raise ValueError("Unsupported encoding; expected UTF-8 or GBK.") from error


def _value(kind: str, token: str) -> Any:
    if kind.lower() == "int":
        return int(token)
    if kind.lower() == "real":
        number = float(token)
        if not math.isfinite(number):
            raise ValueError(f"Non-finite parameter value: {token}")
        return number
    quote = token[0]
    return token[1:-1].replace("\\\\", "\\").replace("\\" + quote, quote)


def _format(kind: str, value: Any, original: str) -> str:
    if kind == "int":
        if isinstance(value, bool):
            raise ValueError("Boolean is not a valid integer parameter.")
        number = int(value)
        if number != value:
            raise ValueError(f"Expected integer value, got {value!r}.")
        return str(number)
    if kind == "real":
        number = float(value)
        if not math.isfinite(number):
            raise ValueError("Real parameter values must be finite.")
        return repr(number)
    quote = original[0]
    return quote + str(value).replace("\\", "\\\\").replace(quote, "\\" + quote) + quote


class ParameterDocument:
    """Keep the original file representation and replace only changed value tokens."""

    def __init__(self, path: Path, raw: bytes, text: str, encoding: str):
        self.path = path
        self.original_bytes = raw
        self.text = text
        self.encoding = encoding
        self.declarations = self._parse(text)

    @classmethod
    def load(cls, path: str | Path) -> "ParameterDocument":
        target = Path(path).expanduser().resolve()
        raw = target.read_bytes()
        text, encoding = _decode(raw)
        return cls(target, raw, text, encoding)

    @staticmethod
    def _parse(text: str) -> tuple[ParameterDeclaration, ...]:
        items = []
        for match in _DECL.finditer(text):
            kind, token = match.group("type"), match.group("value")
            items.append(ParameterDeclaration(
                match.group("name"), kind.lower(), token, _value(kind, token),
                match.start("value"), match.end("value"),
            ))
        return tuple(items)

    def get(self, key: str) -> ParameterDeclaration | None:
        return next((item for item in reversed(self.declarations) if item.name == key), None)

    def render_updates(self, updates: dict[str, Any]) -> bytes:
        """Render validated token replacements without touching the source file."""
        if not updates:
            return self.original_bytes
        replacements = []
        for key, value in updates.items():
            item = self.get(key)
            if item is None:
                raise KeyError(f"Parameter {key!r} is not declared.")
            token = _format(item.declared_type, value, item.raw_value)
            if token != item.raw_value:
                replacements.append((item.value_start, item.value_end, token))
        if not replacements:
            return self.original_bytes
        text = self.text
        for start, end, token in sorted(replacements, reverse=True):
            text = text[:start] + token + text[end:]
        return text.encode(self.encoding)

    def mark_saved(self, raw: bytes) -> None:
        """Update the in-memory snapshot after the corresponding bytes are durable."""
        text = raw.decode(self.encoding)
        self.original_bytes, self.text = raw, text
        self.declarations = self._parse(text)

    def save(self, updates: dict[str, Any]) -> bool:
        if not updates:
            return False
        if self.path.read_bytes() != self.original_bytes:
            raise RuntimeError(
                f"{self.path.name} changed outside the Workbench; reload before saving."
            )
        raw = self.render_updates(updates)
        if raw == self.original_bytes:
            return False
        self.path.write_bytes(raw)
        self.mark_saved(raw)
        return True


def _choice(label: str, value: Any) -> ParameterChoice:
    return ParameterChoice(label, value)


_PARAMETER_SPECS = {
    "simu.txt": (
        ParameterSpec("simutask", "Task type", "OneFLOW task selector.", "choice", "Task",
                      choices=(_choice("Grid generation / conversion", "Grid"), _choice("CFD solve", "Solve"))),
    ),
    "cfd.txt": (
        ParameterSpec("maxSteps", "Maximum iterations", "Maximum solver iteration count.",
                      "integer", "Run control", minimum=1, maximum=2147483647),
        ParameterSpec("maxTime", "Maximum physical time", "Physical-time stopping limit.",
                      "real", "Run control", minimum=0.0, maximum=1.0e12, decimals=8),
        ParameterSpec("time_integral", "Time integration", "Time-integration method.",
                      "choice", "Numerical methods",
                      choices=(_choice("Runge-Kutta (RK)", 1), _choice("LU-SGS", 2))),
        ParameterSpec("ischeme", "Convective flux scheme", "Convective flux discretization.",
                      "choice", "Numerical methods", choices=tuple(
                          _choice(label, value) for label, value in (
                              ("Roe", 1), ("Van Leer", 2), ("Steger-Warming", 3),
                              ("HLLE", 4), ("Lax", 5), ("AUSMP", 6), ("AUSMP+", 7),
                              ("AUSMDV", 8), ("AUSMW", 9), ("AUSMPW", 10), ("SLAU2", 12)
                          )
                      )),
    ),
    "vismodel.txt": (
        ParameterSpec("vismodel", "Viscous model", "Primary viscous-model family.",
                      "choice", "Physical models", choices=tuple(
                          _choice(label, value) for label, value in (
                              ("Euler", 0), ("Navier-Stokes", 1), ("Zero-equation", 2),
                              ("One-equation", 3), ("Two-equation", 4)
                          )
                      )),
        ParameterSpec("visname", "Turbulence / viscosity model", "Named model implementation.",
                      "choice", "Physical models", choices=tuple(
                          _choice(value, value) for value in (
                              "0eq-bl", "1eq-sa", "2eq-kw-menter-sst", "2eq-kw-menter-bsl",
                              "2eq-kw-wilcox-1988", "2eq-kw-wilcox-1998", "2eq-kw-kok-tnt",
                              "2eq-kw-wilcox-2006", "easm-kw-2003", "easm-kw-2005"
                          )
                      )),
    ),
}


def parameter_specs_for_file(path: str | Path) -> tuple[ParameterSpec, ...]:
    """Return metadata backed by checked-in legacy declarations and comments."""
    return _PARAMETER_SPECS.get(Path(path).name.casefold(), ())
