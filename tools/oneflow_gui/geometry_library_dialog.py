""""Browse the format-independent Geometry Asset registry."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (
    QDialog, QLabel, QListWidget, QMessageBox, QPlainTextEdit,
    QPushButton, QSplitter, QVBoxLayout, QWidget,
)


class GeometryLibraryDialog(QDialog):
    """Display registered geometry definitions and their representations."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = Path(repository_root).resolve()
        self.setWindowTitle("Geometry Library")
        self.resize(900, 620)
        root = QVBoxLayout(self)
        heading = QLabel("Geometry Library")
        heading.setStyleSheet("font-size: 22px; font-weight: 650;")
        root.addWidget(heading)
        intro = QLabel(
            "Browse reusable geometry assets independently of file format. "
            "An asset may have analytic or parametric definitions, CAD models, "
            "surface meshes, volume meshes, point clouds or implicit representations."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)
        splitter = QSplitter(self)
        self.asset_list = QListWidget(splitter)
        self.details = QPlainTextEdit(splitter)
        self.details.setReadOnly(True)
        splitter.addWidget(self.asset_list)
        splitter.addWidget(self.details)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 3)
        root.addWidget(splitter, 1)
        self.status = QLabel()
        root.addWidget(self.status)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        root.addWidget(close_button)
        self.asset_list.currentRowChanged.connect(self._show_asset)
        self._assets = self._load_assets()
        for asset in self._assets:
            self.asset_list.addItem(f"{asset['name']}  [{asset['id']}]")
        if self._assets:
            self.asset_list.setCurrentRow(0)
        else:
            self.details.setPlainText("No geometry assets are registered yet.")

    def _load_assets(self) -> list[dict]:
        catalog_path = self._repository_root / "grid" / "geometry" / "catalog.json"
        try:
            catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
            if catalog.get("schema") != "oneflow.geometry-catalog":
                raise ValueError("Unrecognized geometry catalog schema.")
            assets = catalog.get("assets")
            if not isinstance(assets, list):
                raise ValueError("Geometry catalog assets must be a list.")
            return assets
        except (OSError, ValueError, TypeError) as error:
            self.status.setText(f"Geometry catalog unavailable: {error}")
            return []

    def _show_asset(self, row: int) -> None:
        if row < 0 or row >= len(self._assets):
            self.details.clear()
            return
        asset = self._assets[row]
        definition = asset.get("definition", {})
        provenance = asset.get("provenance", {})
        lines = [
            asset.get("name", "Unnamed asset"),
            "=" * len(asset.get("name", "Unnamed asset")),
            f"Asset ID: {asset.get('id', '')}",
            f"Dimension: {asset.get('dimension', 'unknown')}D",
            f"Description: {asset.get('description', '')}",
            "",
            "Definition",
            f"  Kind: {definition.get('kind', 'unknown')}",
        ]
        for equation in definition.get("equations", []):
            lines.append(f"  Equation: {equation}")
        parameters = definition.get("parameters", {})
        if parameters:
            lines.append("  Parameters:")
            lines.extend(f"    {key}: {value}" for key, value in parameters.items())
        if definition.get("notes"):
            lines.append(f"  Notes: {definition['notes']}")
        lines.extend(["", "Provenance", f"  Kind: {provenance.get('kind', 'unknown')}"])
        for key in ("source_case", "tutorial", "source_revision", "notes"):
            if provenance.get(key):
                lines.append(f"  {key.replace('_', ' ').title()}: {provenance[key]}")
        generator = provenance.get("generator", {})
        if generator:
            lines.append(f"  Generator: {generator.get('name', '')}")
            lines.append(f"  Generator source: {generator.get('source', '')}")
            if generator.get("algorithm"):
                lines.append(f"  Algorithm: {generator['algorithm']}")
        lines.extend(["", "Representations"])
        for item in asset.get("representations", []):
            lines.append(
                f"  - {item.get('role', 'unknown')} | {item.get('format', 'unknown')}"
                f" | {item.get('status', 'unspecified')}"
            )
            if item.get("path"):
                lines.append(f"      Path: {item['path']}")
            if item.get("parent_representation"):
                lines.append(f"      Parent: {item['parent_representation']}")
            if item.get("notes"):
                lines.append(f"      Notes: {item['notes']}")
        tutorials = asset.get("tutorials", [])
        if tutorials:
            lines.extend(["", "Tutorials / examples"])
            lines.extend(f"  - {item}" for item in tutorials)
        self.details.setPlainText("\\n".join(lines))
        self.status.setText(
            f"{len(self._assets)} registered assets · metadata only; format support depends on available adapters."
        )
