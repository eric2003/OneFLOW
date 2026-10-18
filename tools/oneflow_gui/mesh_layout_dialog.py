"""Interactive Qt editor for simple structured mesh layout specifications."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from PySide6.QtCore import QByteArray
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from tools.mesh_layout.editor_document import (
    compile_specification,
    load_specification,
    save_specification,
)


class MeshLayoutDialog(QDialog):
    """Edit template geometry and spacing parameters, preview, and compile layouts."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self._spec: dict[str, Any] = {}
        self._field_widgets: dict[tuple[str | int, ...], QLineEdit] = {}

        self.setWindowTitle("Mesh Layout Editor")
        self.resize(1100, 760)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Choose a template, edit geometry and line dimensions, preview the model, "
            "then compile it to the existing OneFLOW layout format."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        template_row = QHBoxLayout()
        template_row.addWidget(QLabel("Template"))
        self.template_combo = QComboBox()
        self.template_combo.addItem("Single rectangle", "rectangle")
        self.template_combo.addItem("Two-block flat plate", "flat_plate")
        self.template_combo.currentIndexChanged.connect(self._load_selected_template)
        template_row.addWidget(self.template_combo, 1)
        root.addLayout(template_row)

        spec_row = QHBoxLayout()
        spec_row.addWidget(QLabel("Specification"))
        self.spec_path_edit = QLineEdit()
        spec_row.addWidget(self.spec_path_edit, 1)
        load_button = QPushButton("Load JSON...")
        load_button.clicked.connect(self._choose_specification)
        spec_row.addWidget(load_button)
        root.addLayout(spec_row)

        body = QHBoxLayout()
        self.fields_area = QScrollArea()
        self.fields_area.setWidgetResizable(True)
        self.fields_panel = QWidget()
        self.fields_form = QFormLayout(self.fields_panel)
        self.fields_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self.fields_area.setWidget(self.fields_panel)
        body.addWidget(self.fields_area, 1)

        self.preview = QSvgWidget()
        self.preview.setMinimumSize(500, 360)
        body.addWidget(self.preview, 2)
        root.addLayout(body, 1)

        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Layout output"))
        self.output_path_edit = QLineEdit()
        output_row.addWidget(self.output_path_edit, 1)
        output_browse = QPushButton("Browse...")
        output_browse.clicked.connect(self._choose_output)
        output_row.addWidget(output_browse)
        root.addLayout(output_row)

        actions = QHBoxLayout()
        save_button = QPushButton("Save specification as...")
        save_button.clicked.connect(self._save_specification_as)
        actions.addWidget(save_button)
        actions.addStretch(1)
        preview_button = QPushButton("Update preview")
        preview_button.clicked.connect(self._update_preview_from_fields)
        actions.addWidget(preview_button)
        compile_button = QPushButton("Compile layout")
        compile_button.setObjectName("primaryButton")
        compile_button.clicked.connect(self._compile_layout)
        actions.addWidget(compile_button)
        close_button = QPushButton("Close")
        close_button.clicked.connect(self.accept)
        actions.addWidget(close_button)
        root.addLayout(actions)

        self.template_combo.setCurrentIndex(0)
        self._load_selected_template()

    def _template_path(self, kind: str) -> Path:
        if kind == "rectangle":
            return self._repository_root / "grid" / "minimal_rectangle" / "rectangle.json"
        return self._repository_root / "grid" / "turbulent_flat_plate" / "flat_plate.json"

    def _default_output_path(self, kind: str) -> Path:
        if kind == "rectangle":
            return self._repository_root / "grid" / "minimal_rectangle" / "grid" / "minimalRectangleLayout.txt"
        return self._repository_root / "grid" / "turbulent_flat_plate" / "grid" / "turbulentFlatPlateLayout.txt"

    def _load_selected_template(self) -> None:
        kind = self.template_combo.currentData()
        path = self._template_path(kind)
        try:
            spec = load_specification(path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Cannot load template", str(error))
            return
        self._spec = spec
        self.spec_path_edit.setText(str(path))
        self.output_path_edit.setText(str(self._default_output_path(kind)))
        self._rebuild_fields()
        self._render_preview()

    def _choose_specification(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load mesh layout specification",
            str(self._repository_root / "grid"),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            spec = load_specification(Path(path))
        except (OSError, ValueError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Invalid specification", str(error))
            return
        kind = spec["type"]
        index = self.template_combo.findData(kind)
        if index < 0:
            QMessageBox.warning(self, "Unsupported template", f"Unsupported type: {kind}")
            return
        self.template_combo.blockSignals(True)
        self.template_combo.setCurrentIndex(index)
        self.template_combo.blockSignals(False)
        self._spec = spec
        self.spec_path_edit.setText(str(Path(path).resolve()))
        self.output_path_edit.setText(str(Path(path).with_suffix(".layout.txt")))
        self._rebuild_fields()
        self._render_preview()

    def _choose_output(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Compile layout to",
            self.output_path_edit.text(),
            "OneFLOW layout (*.txt);;All files (*)",
        )
        if path:
            self.output_path_edit.setText(path)

    def _rebuild_fields(self) -> None:
        while self.fields_form.rowCount():
            self.fields_form.removeRow(0)
        self._field_widgets.clear()

        if self._spec["type"] == "rectangle":
            self._spec.setdefault(
                "boundaries",
                {f"line{line_id}": line_id for line_id in range(1, 5)},
            )
            fields = []
            for point_index in range(4):
                for axis, axis_name in enumerate(("X", "Y", "Z")):
                    fields.append((
                        ("corners", point_index, axis),
                        f"Point {point_index + 1} {axis_name}",
                    ))
            for edge_index in range(4):
                fields.append((("edgePointCounts", edge_index), f"Line {edge_index + 1} point count"))
            fields.append((("spacing",), "Spacing type"))
            for line_id in range(1, 5):
                fields.append((("boundaries", f"line{line_id}"), f"Boundary: Line {line_id}"))
        else:
            fields = [
                (("upstreamDistance",), "Upstream distance"),
                (("plateLength",), "Plate length"),
                (("domainHeight",), "Domain height"),
                (("upstreamPointCount",), "Upstream line point count"),
                (("platePointCount",), "Plate line point count"),
                (("normalPointCount",), "Normal line point count"),
            ]
            for distribution in ("upstream", "plate", "vertical", "topUpstream", "topPlate"):
                for endpoint, title in (("startSpacing", "start spacing"), ("endSpacing", "end spacing")):
                    fields.append((("distributions", distribution, endpoint), f"{distribution} {title}"))
            for boundary in (
                "upstreamBottom", "plate", "left", "interface", "right",
                "topUpstream", "topPlate",
            ):
                fields.append((("boundaries", boundary), f"Boundary: {boundary}"))

        for path, label in fields:
            value = self._get_path(self._spec, path)
            widget = QLineEdit(str(value))
            widget.setToolTip("Use 'auto' for an automatically selected endpoint spacing.")
            self._field_widgets[path] = widget
            self.fields_form.addRow(label, widget)

    @staticmethod
    def _get_path(data: dict[str, Any], path: tuple[str | int, ...]) -> Any:
        value: Any = data
        for key in path:
            value = value[key]
        return value

    @staticmethod
    def _set_path(data: dict[str, Any], path: tuple[str | int, ...], value: Any) -> None:
        target = data
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value

    def _read_fields(self) -> dict[str, Any]:
        spec = json.loads(json.dumps(self._spec))
        for path, widget in self._field_widgets.items():
            raw = widget.text().strip()
            original = self._get_path(spec, path)
            if path[0] == "corners":
                value = float(raw)
            elif path[0] == "edgePointCounts" or path[-1].endswith("PointCount"):
                value = int(raw)
            elif path[0] == "distributions":
                value = "auto" if raw.lower() == "auto" else float(raw)
            elif path[0] in ("plateLength", "upstreamDistance", "domainHeight"):
                value = float(raw)
            elif path == ("spacing",):
                value = raw
            elif path[0] == "boundaries":
                value = raw
            else:
                value = raw if isinstance(original, str) else float(raw)
            self._set_path(spec, path, value)
        return spec

    def _render_preview(self) -> None:
        try:
            _, svg = compile_specification(self._spec)
            self.preview.load(QByteArray(svg.encode("utf-8")))
        except (ValueError, TypeError, KeyError) as error:
            self.preview.load(QByteArray())
            QMessageBox.warning(self, "Cannot preview layout", str(error))

    def _update_preview_from_fields(self) -> None:
        try:
            self._spec = self._read_fields()
            _, svg = compile_specification(self._spec)
            self.preview.load(QByteArray(svg.encode("utf-8")))
        except (ValueError, TypeError, KeyError) as error:
            QMessageBox.warning(self, "Invalid layout parameters", str(error))

    def _save_specification_as(self) -> None:
        try:
            spec = self._read_fields()
            compile_specification(spec)
        except (ValueError, TypeError, KeyError) as error:
            QMessageBox.warning(self, "Invalid layout parameters", str(error))
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save specification as",
            self.spec_path_edit.text(),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            save_specification(Path(path), spec)
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.critical(self, "Cannot save specification", str(error))
            return
        self._spec = spec
        self.spec_path_edit.setText(str(Path(path).resolve()))
        self._rebuild_fields()
        self._render_preview()

    def _compile_layout(self) -> None:
        try:
            spec = self._read_fields()
            layout, svg = compile_specification(spec)
            output_path = Path(self.output_path_edit.text().strip()).expanduser()
            if not str(self.output_path_edit.text()).strip():
                raise ValueError("Choose a layout output path.")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(layout, encoding="utf-8", newline="\n")
        except (OSError, ValueError, TypeError, KeyError) as error:
            QMessageBox.critical(self, "Cannot compile layout", str(error))
            return
        self._spec = spec
        self.preview.load(QByteArray(svg.encode("utf-8")))
        QMessageBox.information(
            self,
            "Layout compiled",
            f"Compiled layout written to:\n{output_path.resolve()}",
        )


__all__ = ["MeshLayoutDialog"]
