"""Persistent settings shared by the OneFLOW desktop tools."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .user_data_paths import UserDataPaths


DEFAULT_REPOSITORY_ROOT = Path("D:/github/OneFLOW")


def workflow_editor_mode_for_run_mode(run_mode: str) -> str | None:
    """Map a saved run preference to the matching workflow editor draft mode."""
    return {
        "multi_case": "multi_directory",
        "dependency_workflow": "case_graph",
    }.get(run_mode)


@dataclass(frozen=True, slots=True)
class ApplicationSettings:
    """User-specific settings that are independent of run workspaces."""

    repository_root: Path = DEFAULT_REPOSITORY_ROOT
    build_configuration: str = "Debug"
    executable_path: str = ""
    run_mode: str = "single_case"

    @classmethod
    def load(cls, settings_file: Path | None = None) -> "ApplicationSettings":
        path = settings_file or UserDataPaths.discover().settings_file
        if not path.is_file():
            return cls()
        with path.open("r", encoding="utf-8") as stream:
            data = json.load(stream)
        repository = data.get("repository_root")
        if not isinstance(repository, str) or not repository.strip():
            return cls()
        configuration = data.get("build_configuration", "Debug")
        if configuration not in {"Debug", "Release"}:
            configuration = "Debug"
        executable = data.get("executable_path", "")
        if not isinstance(executable, str):
            executable = ""
        run_mode = data.get("run_mode", "single_case")
        if run_mode not in {"single_case", "multi_case", "dependency_workflow"}:
            run_mode = "single_case"
        return cls(Path(repository).expanduser(), configuration, executable, run_mode)

    def save(self, settings_file: Path | None = None) -> None:
        path = settings_file or UserDataPaths.discover().settings_file
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "repository_root": str(self.repository_root.expanduser().resolve()),
            "build_configuration": self.build_configuration,
            "executable_path": self.executable_path,
            "run_mode": self.run_mode,
        }
        with path.open("w", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, indent=2)
            stream.write("\n")


def resolve_repository_root(
    settings: ApplicationSettings,
    discovered_root: Path,
) -> Path:
    """Use the configured checkout when present, with a portable launch fallback."""
    configured = settings.repository_root.expanduser()
    if configured.is_dir():
        return configured.resolve()
    return discovered_root.resolve()
