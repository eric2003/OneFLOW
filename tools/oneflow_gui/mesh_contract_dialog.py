"""Mesh contract workbench: inspect a neutral snapshot and its evidence layers."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QPlainTextEdit,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)


def _ensure_tools_on_path(repository_root: Path) -> None:
    tools_root = str((repository_root / "tools").resolve())
    if tools_root not in sys.path:
        sys.path.insert(0, tools_root)


class MeshContractDialog(QDialog):
    """Browse mesh design information and inspect a MeshSnapshot JSON artifact.

    Snapshot validation and converter execution are deliberately reported as
    separate evidence categories. Loading a JSON file never runs a converter.
    """

    _DOMAIN_TEXT = {
        "Overview": (
            "Mesh is the format-neutral geometry and topology domain. It owns "
            "mesh identity, coordinates, ordered connectivity, faces, adjacency, "
            "physical groups and provenance. Format readers are adapters; they "
            "must not make Gmsh-specific types part of Mesh Core."
        ),
        "Principle": (
            "Preserve authoritative data before deriving caches. Node and cell "
            "IDs are identities, not array offsets. Connectivity order carries "
            "orientation information. Physical group names and dimensions are "
            "preserved first; mapping them to solver boundary conditions is an "
            "explicit semantic operation."
        ),
        "Reference": (
            "Primary references: Gmsh MSH 4.x file format specification; CGNS "
            "SIDS and CGNS Mid-Level Library documentation; SU2 mesh format; "
            "OneFLOW mesh architecture documentation. External format manuals "
            "are references, not proof of adapter conformance."
        ),
        "Architecture": (
            "Format reader -> neutral MeshSnapshot -> contract validation -> "
            "explicit semantic mapping -> production C++ Mesh Core -> writer. "
            "The Python validator is a tooling-side contract checker only. "
            "Production C++ snapshot remains the target and JSON is not a runtime "
            "dependency."
        ),
        "Implementation": (
            "MeshSnapshot v1 currently describes dimensioned nodes, cells, faces "
            "with owner/neighbour IDs, dimensioned groups and source provenance. "
            "The validator checks identity, finite coordinates, element arity "
            "and dimension, connectivity references, face ownership and group "
            "entity kinds. It does not infer orientation or boundary physics."
        ),
        "Code": (
            "Contract validator: tools/mesh_contracts/snapshot.py\n"
            "Contract tests: tools/mesh_contracts/test_snapshot.py\n"
            "Gmsh adapter work is tracked separately; do not treat the current "
            "placeholder output as a production CGNS writer."
        ),
        "Documentation": (
            "doc/source/mesh_architecture.rst\n"
            "tools/mesh_contracts/README.md\n"
            "Gmsh format/interface notes are maintained with the geometry "
            "conversion documentation."
        ),
    }

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        _ensure_tools_on_path(self._repository_root)
        from mesh_contracts.snapshot import MeshSnapshotError, validate_snapshot

        self._snapshot_error = MeshSnapshotError
        self._validate_snapshot = validate_snapshot
        self._snapshot: dict[str, Any] | None = None
        self.setWindowTitle("Contract Workbench - Mesh")
        self.resize(1180, 780)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Mesh has its own topology and interchange knowledge model. "
            "Snapshot validation is not converter execution."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        actions = QHBoxLayout()
        self.open_snapshot_button = QPushButton("Open MeshSnapshot JSON...")
        self.open_snapshot_button.clicked.connect(self._open_snapshot)
        actions.addWidget(self.open_snapshot_button)
        self.snapshot_path_label = QLabel("No snapshot loaded")
        self.snapshot_path_label.setWordWrap(True)
        actions.addWidget(self.snapshot_path_label, 1)
        root.addLayout(actions)

        self.tabs = QTabWidget()
        for name, description in self._DOMAIN_TEXT.items():
            if name == "Code":
                description = self._read_repository_files(
                    ("tools/mesh_contracts/snapshot.py", "codes/geometry/src/GridConversion.cpp"),
                    max_chars=24000,
                )
            elif name == "Documentation":
                description = self._read_repository_files(
                    ("doc/source/mesh_architecture.rst", "doc/source/gmsh_interface.rst")
                )
            view = QPlainTextEdit(description)
            view.setReadOnly(True)
            view.setLineWrapMode(
                QPlainTextEdit.LineWrapMode.NoWrap
                if name == "Code"
                else QPlainTextEdit.LineWrapMode.WidgetWidth
            )
            self.tabs.addTab(view, name)

        topology_page = QWidget()
        self.topology_page = topology_page
        topology_layout = QVBoxLayout(topology_page)
        self.entity_selector = QComboBox()
        self.entity_selector.addItems(["Nodes", "Cells", "Faces", "Groups"])
        self.entity_selector.currentTextChanged.connect(self._show_entity_table)
        topology_layout.addWidget(self.entity_selector)
        self.topology_summary = QLabel("Load a MeshSnapshot to inspect topology.")
        topology_layout.addWidget(self.topology_summary)
        self.topology_table = QTableWidget()
        self.topology_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.topology_table.setAlternatingRowColors(True)
        topology_layout.addWidget(self.topology_table, 1)
        self.tabs.addTab(topology_page, "Topology")

        evidence_page = QWidget()
        evidence_layout = QVBoxLayout(evidence_page)
        evidence_layout.addWidget(QLabel("MeshSnapshot contract validation"))
        self.validation_view = QPlainTextEdit()
        self.validation_view.setReadOnly(True)
        self.validation_view.setPlainText("NOT RUN — no snapshot loaded.")
        evidence_layout.addWidget(self.validation_view, 1)
        evidence_layout.addWidget(QLabel("Format conversion execution"))
        self.conversion_view = QPlainTextEdit()
        self.conversion_view.setReadOnly(True)
        self.conversion_view.setPlainText(
            "NOT RUN — loading or validating a snapshot does not execute a "
            "Gmsh/CGNS/OneFLOW converter. A future converter run will record "
            "its command, exit status, input/output fingerprints and diagnostics here."
        )
        evidence_layout.addWidget(self.conversion_view, 1)
        self.tabs.addTab(evidence_page, "Verification / Evidence")

        root.addWidget(self.tabs, 1)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        root.addWidget(close_button)

    def _read_repository_files(
        self, paths: tuple[str, ...], *, max_chars: int | None = None
    ) -> str:
        """Read declared source/documentation files without allowing root escape."""
        root = self._repository_root.resolve()
        sections: list[str] = []
        for relative in paths:
            sections.extend((f"File: {relative}", "=" * (6 + len(relative)), ""))
            candidate = (root / relative).resolve()
            try:
                candidate.relative_to(root)
                content = candidate.read_text(encoding="utf-8")
            except (ValueError, OSError, UnicodeError) as error:
                sections.append(f"Unavailable: {type(error).__name__}: {error}")
                sections.append("")
                continue
            if max_chars is not None and len(content) > max_chars:
                content = content[:max_chars] + "\n\n[Preview truncated; open the source file for the complete content.]"
            sections.append(content.rstrip())
            sections.append("")
        return "\n".join(sections).rstrip()

    def _open_snapshot(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open MeshSnapshot",
            str(self._repository_root),
            "MeshSnapshot JSON (*.json);;All files (*)",
        )
        if path:
            self._load_snapshot_file(Path(path))

    def _load_snapshot_file(self, path: Path) -> bool:
        """Load and validate one manifest, keeping validation separate from conversion."""
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            self._validate_snapshot(payload)
        except (OSError, UnicodeError, json.JSONDecodeError, self._snapshot_error) as exc:
            self._snapshot = None
            self.snapshot_path_label.setText(str(path))
            self.validation_view.setPlainText(f"FAILED — {type(exc).__name__}: {exc}")
            self.topology_summary.setText("Snapshot is unavailable because validation failed.")
            self.topology_table.setRowCount(0)
            self.topology_table.setColumnCount(0)
            self.tabs.setCurrentWidget(self.tabs.widget(self.tabs.count() - 1))
            return False
        self._snapshot = payload
        self.snapshot_path_label.setText(str(path))
        self.validation_view.setPlainText(
            "PASSED — MeshSnapshot schema and declared invariants validated.\n"
            f"Dimension: {payload['dimension']}\n"
            f"Nodes: {len(payload['nodes'])}\n"
            f"Cells: {len(payload['cells'])}\n"
            f"Faces: {len(payload['faces'])}\n"
            f"Groups: {len(payload['groups'])}\n"
            f"Source: {payload['provenance']['sourceFormat']}\n\n"
            "This result validates the manifest only; it does not prove geometric "
            "quality, face orientation, complete incidence, boundary-condition "
            "semantics, or converter correctness."
        )
        self._show_entity_table(self.entity_selector.currentText())
        self.tabs.setCurrentWidget(self.topology_page)
        return True

    def _show_entity_table(self, entity: str) -> None:
        if self._snapshot is None:
            return
        key = {"Nodes": "nodes", "Cells": "cells", "Faces": "faces", "Groups": "groups"}[entity]
        records = self._snapshot[key]
        columns = {
            "Nodes": ("id", "coordinates"),
            "Cells": ("id", "type", "nodeIds"),
            "Faces": ("id", "nodeIds", "ownerCellId", "neighbourCellId"),
            "Groups": ("id", "name", "kind", "dimension", "memberKind", "memberIds"),
        }[entity]
        self.topology_table.clear()
        self.topology_table.setColumnCount(len(columns))
        self.topology_table.setHorizontalHeaderLabels(list(columns))
        self.topology_table.setRowCount(len(records))
        for row, record in enumerate(records):
            for column, name in enumerate(columns):
                value = record.get(name)
                if isinstance(value, (dict, list)):
                    value = json.dumps(value, ensure_ascii=True)
                self.topology_table.setItem(row, column, QTableWidgetItem(str(value)))
        self.topology_table.resizeColumnsToContents()
        self.topology_summary.setText(
            f"{entity}: {len(records)} record(s). Source IDs and connectivity order are shown as stored."
        )
