"""Persistent settings shared by the OneFLOW desktop tools."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from .user_data_paths import UserDataPaths


DEFAULT_REPOSITORY_ROOT = Path("D:/github/OneFLOW")


@dataclass(frozen=True, slots=True)
class ApplicationSettings:
    """User-specific settings that are independent of run workspaces."""

    repository_root: Path = DEFAULT_REPOSITORY_ROOT
    build_configuration: str = "Debug"
    executable_path: str = ""

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
        return cls(Path(repository).expanduser(), configuration, executable)

    def save(self, settings_file: Path | None = None) -> None:
        path = settings_file or UserDataPaths.discover().settings_file
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "repository_root": str(self.repository_root.expanduser().resolve()),
            "build_configuration": self.build_configuration,
            "executable_path": self.executable_path,
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
