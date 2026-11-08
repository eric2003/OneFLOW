"""Query and display dependency information from the running OneFLOW binary."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

from PySide6.QtWidgets import (
    QDialog,
    QHeaderView,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)


def query_runtime_environment(
    executable: str | Path,
    *,
    timeout: float = 15.0,
) -> dict[str, Any]:
    """Ask the selected OneFLOW executable for its runtime environment JSON."""
    path = Path(executable).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"OneFLOW executable does not exist: {path}")

    command = [str(path), "--environment", "--format", "json"]
    try:
        result = subprocess.run(
            command,
            cwd=path.parent,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"OneFLOW environment query exceeded {timeout:g} seconds."
        ) from error

    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip()
        raise RuntimeError(
            f"OneFLOW environment query failed with exit code "
            f"{result.returncode}: {detail}"
        )

    try:
        document = json.loads(result.stdout)
    except json.JSONDecodeError as error:
        raise ValueError(
            "OneFLOW returned invalid environment JSON."
        ) from error

    if not isinstance(document, dict) or document.get("schema_version") != 1:
        raise ValueError("Unsupported OneFLOW environment information schema.")
    if not isinstance(document.get("oneflow"), dict):
        raise ValueError("OneFLOW environment JSON has no oneflow object.")
    if not isinstance(document.get("build"), dict):
        raise ValueError("OneFLOW environment JSON has no build object.")
    if not isinstance(document.get("dependencies"), list):
        raise ValueError("OneFLOW environment JSON has no dependencies list.")

    return document


class RuntimeEnvironmentDialog(QDialog):
    """Display runtime versions and the modules loaded by a OneFLOW process."""

    def __init__(
        self,
        executable: str | Path,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._executable = Path(executable).expanduser().resolve()

        self.setWindowTitle("OneFLOW Runtime Environment")
        self.resize(1120, 520)

        layout = QVBoxLayout(self)
        self.summary_label = QLabel("Querying the selected OneFLOW executable...")
        self.summary_label.setWordWrap(True)
        layout.addWidget(self.summary_label)

        self.table = QTableWidget(0, 6, self)
        self.table.setHorizontalHeaderLabels(
            [
                "Dependency",
                "Runtime version",
                "Header version",
                "Version consistency",
                "Status",
                "Loaded module",
            ]
        )
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectRows
        )
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.ResizeToContents
        )
        self.table.horizontalHeader().setStretchLastSection(True)
        layout.addWidget(self.table, 1)

        button_layout = QHBoxLayout()
        button_layout.addStretch(1)
        refresh_button = QPushButton("Refresh")
        refresh_button.clicked.connect(self.refresh)
        button_layout.addWidget(refresh_button)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        button_layout.addWidget(close_button)
        layout.addLayout(button_layout)

        self.refresh()

    def refresh(self) -> None:
        """Refresh the table from the executable selected in Workbench settings."""
        try:
            document = query_runtime_environment(self._executable)
        except (OSError, RuntimeError, ValueError) as error:
            self.summary_label.setText(
                f"Could not query runtime environment: {error}"
            )
            self.table.setRowCount(0)
            return

        oneflow = document["oneflow"]
        build = document["build"]
        self.summary_label.setText(
            f"OneFLOW {oneflow.get('version', 'unknown')} | "
            f"{build.get('compiler', 'unknown')} "
            f"{build.get('compiler_version', '')} | "
            f"{build.get('configuration', 'unknown')} | "
            f"Executable: {self._executable}"
        )

        dependencies = document["dependencies"]
        self.table.setRowCount(len(dependencies))
        for row, dependency in enumerate(dependencies):
            runtime = dependency.get("runtime_version") or {}
            header = dependency.get("header_version") or {}
            runtime_value = runtime.get("value")
            runtime_source = runtime.get("source")
            runtime_text = (
                f"{runtime_value} ({runtime_source})"
                if runtime_value
                else "Unavailable"
            )
            header_text = header.get("value") or "—"
            compatibility = dependency.get("version_compatibility") or {}
            compatibility_status = compatibility.get("status", "unknown")
            compatibility_basis = compatibility.get("basis", "")
            compatibility_text = (
                f"{compatibility_status} ({compatibility_basis})"
                if compatibility_basis
                else compatibility_status
            )
            compatibility_detail = compatibility.get("detail", "")
            runtime_library = dependency.get("runtime_library") or {}
            module_path = runtime_library.get("module_path") or "—"
            linkage = runtime_library.get("linkage")
            if linkage:
                module_text = f"{linkage}: {module_path}"
            else:
                module_text = module_path

            values = (
                dependency.get("name", "Unknown"),
                runtime_text,
                header_text,
                compatibility_text,
                dependency.get("status", "unknown"),
                module_text,
            )
            for column, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                if column == 4:
                    item.setToolTip(str(compatibility_detail))
                elif column == 5:
                    item.setToolTip(str(module_path))
                self.table.setItem(row, column, item)

        self.table.resizeRowsToContents()
