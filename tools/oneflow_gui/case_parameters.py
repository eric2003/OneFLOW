"""Lossless reader for legacy OneFLOW parameter files."""
from __future__ import annotations
import re
from dataclasses import dataclass
from pathlib import Path

_DECL = re.compile(r"(?m)^(?P<type>int|Real|string)\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)\s*=\s*(?P<value>[^;\r\n]*?)(?P<tail>\s*(?:;\s*)?(?://.*)?)$")

@dataclass(frozen=True, slots=True)
class ParameterDeclaration:
    name: str
    declared_type: str
    raw_value: str
    line: int
    start: int
    end: int

@dataclass(frozen=True, slots=True)
class ParameterDocument:
    path: Path
    encoding: str
    text: str
    declarations: tuple[ParameterDeclaration, ...]

    @classmethod
    def read(cls, path: str | Path) -> "ParameterDocument":
        path = Path(path)
        data = path.read_bytes()
        encoding = "utf-8-sig" if data.startswith(b"\xef\xbb\xbf") else "utf-8"
        if encoding == "utf-8":
            try:
                data.decode("utf-8")
            except UnicodeDecodeError:
                encoding = "gb18030"
        text = data.decode(encoding)
        items = tuple(ParameterDeclaration(m["name"], m["type"], m["value"].strip(),
            text.count("\n", 0, m.start()) + 1, m.start("value"), m.end("value"))
            for m in _DECL.finditer(text))
        return cls(path, encoding, text, items)

    def values(self, name: str) -> tuple[ParameterDeclaration, ...]:
        return tuple(item for item in self.declarations if item.name == name)

    def replace_value(self, name: str, value: str | int | float, occurrence: int = -1) -> bytes:
        items = self.values(name)
        if not items:
            raise KeyError(name)
        try:
            item = items[occurrence]
        except IndexError as error:
            raise IndexError(name) from error
        value = str(value)
        if "\n" in value or "\r" in value:
            raise ValueError("Parameter value must be a single line.")
        return (self.text[:item.start] + value + self.text[item.end:]).encode(self.encoding)

    def save_value(self, name: str, value: str | int | float, occurrence: int = -1) -> None:
        self.path.write_bytes(self.replace_value(name, value, occurrence))

def discover_parameter_documents(control_file: str | Path) -> tuple[ParameterDocument, ...]:
    manifest = Path(control_file)
    documents = []
    for line in manifest.read_text(encoding="utf-8-sig").splitlines():
        entry = line.split("//", 1)[0].strip().rstrip(";").strip().strip('"').strip("'")
        if not entry:
            continue
        candidate = (manifest.parent / entry).resolve()
        try:
            candidate.relative_to(manifest.parent.resolve())
        except ValueError as error:
            raise ValueError(f"Parameter path escapes case directory: {entry}") from error
        if candidate.is_file():
            documents.append(ParameterDocument.read(candidate))
    return tuple(documents)
