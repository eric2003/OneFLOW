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
        return cls(Path(repository).expanduser())

    def save(self, settings_file: Path | None = None) -> None:
        path = settings_file or UserDataPaths.discover().settings_file
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"repository_root": str(self.repository_root.expanduser().resolve())}
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
