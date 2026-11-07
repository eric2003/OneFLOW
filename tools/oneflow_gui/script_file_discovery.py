"""Discover parameter declarations and validate optional filename-only manifests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True, slots=True)
class ScriptFileDiscovery:
    """Files to parse, their discovery mode, and non-fatal manifest diagnostics."""

    files: tuple[Path, ...]
    mode: str
    diagnostics: tuple[str, ...] = ()


def _read_text(path: Path) -> str:
    data = path.read_bytes()
    for encoding in ("utf-8-sig", "gbk"):
        try:
            return data.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise UnicodeError(f"unsupported control-file encoding: {path}")


def _strip_comment(line: str) -> str:
    return re.split(r"#|//", line, maxsplit=1)[0].strip()


def _manifest_name(line: str) -> str | None:
    line = _strip_comment(line)
    if not line:
        return ""
    match = re.fullmatch(r"""(?:"([^"]+)"|'([^']+)'|([^\s;]+))\s*;?""", line)
    if match is None:
        return None
    name = next((part for part in match.groups() if part is not None), "")
    candidate = Path(name)
    if (
        not name
        or candidate.is_absolute()
        or candidate.name != name
        or candidate.suffix != ".txt"
        or name in {".", "..", "control.txt"}
    ):
        return None
    return name


def _has_parameter_declaration(path: Path) -> bool:
    try:
        content = _read_text(path)
    except (OSError, UnicodeError):
        # Keep unreadable files in the ingestion list so the parser surfaces
        # the actual file error instead of silently hiding the file.
        return True
    return any("=" in _strip_comment(line) for line in content.splitlines())


def discover_script_parameter_files(
    script_directory: str | Path,
) -> ScriptFileDiscovery:
    """Discover .txt files and classify control.txt by its actual contents.

    Manifest entries preserve their historical precedence. Existing files
    omitted from a manifest are appended in filename order and reported.
    Missing and duplicate entries are reported but never loaded.
    """
    script = Path(script_directory).expanduser().resolve()
    if not script.is_dir():
        return ScriptFileDiscovery(
            (), "directory", (f"Script directory does not exist: {script}",)
        )

    control = script / "control.txt"
    actual = tuple(sorted(
        (path.resolve() for path in script.glob("*.txt") if path.is_file()),
        key=lambda path: path.name,
    ))
    if not control.is_file():
        parameter_files = tuple(path for path in actual if _has_parameter_declaration(path))
        auxiliary = tuple(path for path in actual if path not in parameter_files)
        diagnostics = (
            "script/control.txt is absent; reading parameter files discovered in script/.",
            *(f"Auxiliary script file {path.name!r} has no parameter declarations and is not loaded as a parameter file." for path in auxiliary),
        )
        return ScriptFileDiscovery(parameter_files, "directory", diagnostics)

    try:
        content = _read_text(control)
    except (OSError, UnicodeError) as error:
        parameter_files = tuple(
            path for path in actual
            if path.name == "control.txt" or _has_parameter_declaration(path)
        )
        auxiliary = tuple(
            path for path in actual
            if path.name != "control.txt" and path not in parameter_files
        )
        diagnostics = (
            f"control.txt could not be classified ({error}); reading it as a parameter file.",
            *(f"Auxiliary script file {path.name!r} has no parameter declarations and is not loaded as a parameter file." for path in auxiliary),
        )
        return ScriptFileDiscovery(parameter_files, "parameter-file", diagnostics)

    names: list[str] = []
    is_manifest = True
    for line in content.splitlines():
        name = _manifest_name(line)
        if name is None:
            is_manifest = False
            break
        if name:
            names.append(name)
    if not names:
        is_manifest = False

    if not is_manifest:
        parameter_files = tuple(
            path for path in actual
            if path.name == "control.txt" or _has_parameter_declaration(path)
        )
        auxiliary = tuple(
            path for path in actual
            if path.name != "control.txt" and path not in parameter_files
        )
        diagnostics = (
            "script/control.txt is not a filename-only manifest; reading it as a parameter file.",
            *(f"Auxiliary script file {path.name!r} has no parameter declarations and is not loaded as a parameter file." for path in auxiliary),
        )
        return ScriptFileDiscovery(parameter_files, "parameter-file", diagnostics)

    all_by_name = {path.name: path for path in actual if path.name != "control.txt"}
    actual_by_name = {
        name: path for name, path in all_by_name.items()
        if _has_parameter_declaration(path)
    }
    diagnostics: list[str] = [
        f"Auxiliary script file {name!r} has no parameter declarations and is not loaded as a parameter file."
        for name, path in all_by_name.items()
        if name not in actual_by_name
    ]
    selected: list[Path] = []
    declared: set[str] = set()
    for name in names:
        if name in declared:
            diagnostics.append(f"script/control.txt lists {name!r} more than once.")
            continue
        declared.add(name)
        path = actual_by_name.get(name)
        if path is None:
            if name in all_by_name:
                diagnostics.append(
                    f"script/control.txt lists {name!r}, but it has no parameter declarations and is not loaded."
                )
            else:
                diagnostics.append(f"script/control.txt lists missing parameter file {name!r}.")
            continue
        selected.append(path)

    for path in actual:
        if (
            path.name == "control.txt"
            or path.name in declared
            or path.name not in actual_by_name
        ):
            continue
        diagnostics.append(
            f"Parameter file {path.name!r} exists in script/ but is not listed in control.txt; it will still be read."
        )
        selected.append(path)

    return ScriptFileDiscovery(tuple(selected), "manifest", tuple(diagnostics))
