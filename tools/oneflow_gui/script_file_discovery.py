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
        or "/" in name
        or "\\" in name
        or ":" in name
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

    Actual parameter files are always loaded in filename order. A filename-only
    manifest is checked for missing, extra, duplicate, and order differences;
    it never determines ingestion or precedence.
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
            f"control.txt could not be classified ({error}); reading it as an ordinary parameter file.",
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
            "script/control.txt is not a filename-only manifest; reading it as an ordinary parameter file.",
            *(f"Auxiliary script file {path.name!r} has no parameter declarations and is not loaded as a parameter file." for path in auxiliary),
        )
        return ScriptFileDiscovery(parameter_files, "parameter-file", diagnostics)

    # Real files define ingestion. A filename-only control.txt is only
    # checked against the directory and never controls which files are read.
    all_by_name = {path.name: path for path in actual if path.name != "control.txt"}
    actual_by_name = {
        name: path for name, path in all_by_name.items()
        if _has_parameter_declaration(path)
    }
    actual_names = tuple(actual_by_name)
    diagnostics: list[str] = [
        f"Auxiliary script file {name!r} has no parameter declarations and is not loaded as a parameter file."
        for name in all_by_name
        if name not in actual_by_name
    ]

    declared: set[str] = set()
    existing_manifest_names: list[str] = []
    for name in names:
        if name in declared:
            diagnostics.append(f"script/control.txt lists {name!r} more than once.")
            continue
        declared.add(name)
        if name not in actual_by_name:
            if name in all_by_name:
                diagnostics.append(
                    f"script/control.txt lists {name!r}, but it has no parameter declarations and is not loaded."
                )
            else:
                diagnostics.append(
                    f"WARNING: script/control.txt lists missing parameter file {name!r}; "
                    "only files present in script/ will be read."
                )
            continue
        existing_manifest_names.append(name)

    for name in actual_names:
        if name not in declared:
            diagnostics.append(
                f"Parameter file {name!r} exists in script/ but is not listed in control.txt; it will still be read."
            )

    if tuple(existing_manifest_names) != actual_names:
        diagnostics.append(
            "Parameter-file order in script/control.txt differs from the actual "
            "script/ filename order; real parameter files will be read in filename order."
        )

    selected = [actual_by_name[name] for name in actual_names]
    return ScriptFileDiscovery(tuple(selected), "manifest", tuple(diagnostics))
