"""Centralized, platform-independent paths for OneFLOW user data."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class UserDataPaths:
    """Resolve user-owned application data without creating directories."""

    root: Path

    @classmethod
    def discover(
        cls,
        *,
        home: Path | None = None,
        environ: dict[str, str] | None = None,
    ) -> "UserDataPaths":
        """Use ONEFLOW_HOME when configured, otherwise the user's home folder."""
        environment = os.environ if environ is None else environ
        configured_root = environment.get("ONEFLOW_HOME", "").strip()
        if configured_root:
            root = Path(configured_root).expanduser()
        else:
            root = (home if home is not None else Path.home()) / ".oneflow"
        return cls(root.resolve())

    @property
    def settings_file(self) -> Path:
        return self.root / "settings.json"

    @property
    def mesh_layout_workspace(self) -> Path:
        """Keep editable mesh specifications and generated layouts outside the repository."""
        return self.root / "mesh-layout-workspace"

    @property
    def mesh_layout_projects(self) -> Path:
        return self.mesh_layout_workspace / "projects"

    @property
    def mesh_layout_specifications(self) -> Path:
        return self.mesh_layout_workspace / "specifications"

    @property
    def mesh_layout_outputs(self) -> Path:
        return self.mesh_layout_workspace / "layouts"

    @property
    def workflow_data(self) -> Path:
        """Keep the existing data location so saved workflows and history survive."""
        return self.root / "workflow-data"

    @property
    def workflow_definitions(self) -> Path:
        return self.workflow_data / "definitions"

    @property
    def workflow_runs(self) -> Path:
        return self.workflow_data / "runs"

    @property
    def executions(self) -> Path:
        return self.workflow_data / "executions"

    @property
    def logs(self) -> Path:
        return self.root / "logs"
