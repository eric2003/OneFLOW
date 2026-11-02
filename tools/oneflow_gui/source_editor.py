"""Explicit source-editor discovery and process launch for Workbench files."""
from __future__ import annotations
import os
import shutil
import subprocess
from pathlib import Path
from typing import Sequence

EDITOR_CHOICES = ("vscode", "vs2026", "vs2022", "notepad++", "system")

def _existing_file(candidates: Sequence[Path]) -> Path | None:
    return next((path.resolve() for path in candidates if path.is_file()), None)

def _visual_studio_executable(editor: str) -> Path | None:
    root = Path(os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)"))
    vswhere = root / "Microsoft Visual Studio" / "Installer" / "vswhere.exe"
    if not vswhere.is_file():
        return None
    version = "[18.0,19.0)" if editor == "vs2026" else "[17.0,18.0)"
    result = subprocess.run([str(vswhere), "-latest", "-products", "*", "-version", version,
        "-find", r"Common7\IDE\devenv.exe"], capture_output=True, text=True, check=False,
        creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
    return _existing_file(tuple(Path(line.strip()) for line in result.stdout.splitlines() if line.strip()))

def resolve_editor_executable(editor: str) -> Path:
    if editor not in EDITOR_CHOICES:
        raise ValueError(f"Unsupported source editor: {editor}")
    if editor == "system":
        raise ValueError("Choose an explicit editor for source files.")
    if editor in ("vs2026", "vs2022"):
        resolved = _visual_studio_executable(editor) if os.name == "nt" else None
    elif editor == "vscode":
        command = shutil.which("code")
        roots = [Path(os.environ[key]) for key in ("LOCALAPPDATA", "ProgramFiles") if os.environ.get(key)]
        candidates = [root / "Programs/Microsoft VS Code/Code.exe" for root in roots[:1]]
        candidates += [root / "Microsoft VS Code/Code.exe" for root in roots[1:]]
        resolved = Path(command).resolve() if command else _existing_file(candidates)
    else:
        command = shutil.which("notepad++") or shutil.which("notepad++.exe")
        roots = [Path(os.environ[key]) for key in ("ProgramFiles", "ProgramFiles(x86)") if os.environ.get(key)]
        candidates = [root / "Notepad++/notepad++.exe" for root in roots]
        resolved = Path(command).resolve() if command else _existing_file(candidates)
    if resolved is None:
        label = {"vscode":"Visual Studio Code", "vs2026":"Visual Studio 2026",
                 "vs2022":"Visual Studio 2022", "notepad++":"Notepad++"}[editor]
        raise FileNotFoundError(f"{label} was not found. Install it or choose another editor in OneFLOW Configuration.")
    return resolved

def open_source_files(editor: str, files: Sequence[Path], *, working_directory: Path) -> subprocess.Popen:
    paths = tuple(Path(path).resolve() for path in files)
    if not paths:
        raise ValueError("No source files were selected.")
    if any(not path.is_file() for path in paths):
        raise FileNotFoundError("One or more selected files no longer exist.")
    executable = resolve_editor_executable(editor)
    arguments = [str(executable)]
    if editor == "vscode":
        arguments.append("--reuse-window")
    arguments.extend(str(path) for path in paths)
    return subprocess.Popen(arguments, cwd=str(Path(working_directory).resolve()))
