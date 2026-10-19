"""Interactive Qt editor for simple structured mesh layout specifications."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from PySide6.QtCore import QByteArray, Qt
from PySide6.QtSvgWidgets import QSvgWidget
from PySide6.QtWidgets import (
    QCheckBox,
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

from .user_data_paths import UserDataPaths
from .mesh_layout_projects import ensure_template_project


BOUNDARY_OPTIONS = (
    ("none", 0, "No physical boundary condition is applied."),
    ("interface", -1, "Connects neighboring mesh blocks; not an external boundary."),
    ("extrapolation", 1, "Extrapolates the boundary state from nearby interior cells."),
    ("wall", 2, "Solid wall; use for the plate or another solid surface."),
    ("symmetry", 3, "Symmetry plane; normal flow is zero."),
    ("farfield", 4, "External boundary representing the far-field flow."),
    ("inflow", 5, "Inlet where the incoming flow state is prescribed."),
    ("outflow", 6, "Outlet where flow leaves the computational domain."),
    ("interior", 100, "Internal region marker; normally not an outer boundary."),
    ("overset", 1000, "Boundary used for overlapping or overset grid connectivity."),
)
BOUNDARY_NAMES_BY_CODE = {code: name for name, code, _ in BOUNDARY_OPTIONS}


class MeshLayoutDialog(QDialog):
    """Edit template geometry and spacing parameters, preview, and compile layouts."""

    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        self._spec: dict[str, Any] = {}
        self._project_root: Path | None = None
        self._field_widgets: dict[tuple[str | int, ...], QLineEdit | QComboBox] = {}

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
        self.template_combo.addItem("Single rectangle with vertical midpoints", "rectangle_vertical_midpoints")
        self.template_combo.addItem("Two-block turbulent flat plate", "flat_plate")
        self.template_combo.addItem("Two-block laminar flat plate", "laminar_flat_plate")
        self.template_combo.addItem("Two-block laminar flat plate (10 edges)", "laminar_flat_plate_10_edges")
        self.template_combo.currentIndexChanged.connect(self._load_selected_template)
        template_row.addWidget(self.template_combo, 1)
        root.addLayout(template_row)

        boundary_display_row = QHBoxLayout()
        self.show_boundary_codes = QCheckBox("Show boundary codes")
        self.show_boundary_codes.setToolTip(
            "Off: show concise names such as Wall. On: show names and legacy codes."
        )
        self.show_boundary_codes.toggled.connect(self._refresh_boundary_labels)
        boundary_display_row.addWidget(self.show_boundary_codes)
        boundary_display_row.addWidget(
            QLabel("Hover over a boundary option to read its explanation.")
        )
        self.reset_boundary_button = QPushButton("Reset flat-plate boundaries")
        self.reset_boundary_button.clicked.connect(
            self._reset_flat_plate_boundary_defaults
        )
        boundary_display_row.addWidget(self.reset_boundary_button)
        boundary_display_row.addStretch(1)
        root.addLayout(boundary_display_row)

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
        save_button = QPushButton("Save specification")
        save_button.clicked.connect(self._save_specification)
        actions.addWidget(save_button)
        save_as_button = QPushButton("Save specification as...")
        save_as_button.clicked.connect(self._save_specification_as)
        actions.addWidget(save_as_button)
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

    def _template_project_source(self, kind: str) -> Path:
        folder = {"rectangle": "minimal_rectangle", "rectangle_vertical_midpoints": "minimal_rectangle_vertical_midpoints", "flat_plate": "turbulent_flat_plate", "laminar_flat_plate": "laminar_flat_plate", "laminar_flat_plate_10_edges": "laminar_flat_plate_10_edges"}[kind]
        return self._repository_root / "grid" / folder

    def _template_specification_name(self, kind: str) -> str:
        return {"rectangle": "rectangle.json", "rectangle_vertical_midpoints": "rectangle_vertical_midpoints.json", "flat_plate": "flat_plate.json", "laminar_flat_plate": "flat_plate.json", "laminar_flat_plate_10_edges": "flat_plate.json"}[kind]

    def _layout_output_name(self, kind: str) -> str:
        return {
            "rectangle": "minimalRectangleLayout.txt",
            "rectangle_vertical_midpoints": "minimalRectangleVerticalMidpointsLayout.txt",
            "laminar_flat_plate": "laminarFlatPlateLayout.txt",
            "laminar_flat_plate_10_edges": "laminarFlatPlate10EdgeLayout.txt",
            "flat_plate": "turbulentFlatPlateLayout.txt",
        }[kind]

    def _workspace_paths(self) -> UserDataPaths:
        return UserDataPaths.discover()

    def _default_specification_path(self, kind: str) -> Path:
        if self._project_root is not None:
            return self._project_root / self._template_specification_name(kind)
        return self._workspace_paths().mesh_layout_projects / (
            {"rectangle": "minimal_rectangle", "rectangle_vertical_midpoints": "minimal_rectangle_vertical_midpoints", "flat_plate": "turbulent_flat_plate", "laminar_flat_plate": "laminar_flat_plate", "laminar_flat_plate_10_edges": "laminar_flat_plate_10_edges"}[kind]
        ) / self._template_specification_name(kind)

    def _default_output_path(self, kind: str, project_root: Path | None = None) -> Path:
        root = project_root or self._project_root
        if root is None:
            root = self._workspace_paths().mesh_layout_projects / (
                {"rectangle": "minimal_rectangle", "rectangle_vertical_midpoints": "minimal_rectangle_vertical_midpoints", "flat_plate": "turbulent_flat_plate",
                 "laminar_flat_plate": "laminar_flat_plate", "laminar_flat_plate_10_edges": "laminar_flat_plate_10_edges"}[kind]
            )
        return root / "grid" / self._layout_output_name(kind)

    def _require_outside_repository(self, path: Path) -> Path:
        resolved_path = path.expanduser().resolve()
        try:
            resolved_path.relative_to(self._repository_root)
        except ValueError:
            return resolved_path
        raise ValueError(
            "Mesh Layout Editor will not write into the OneFLOW source tree. "
            "Open a copied template project or choose a project outside the repository."
        )

    def _load_selected_template(self) -> None:
        kind = self.template_combo.currentData()
        try:
            self._project_root = ensure_template_project(
                self._repository_root,
                self._workspace_paths().mesh_layout_projects,
                kind,
            )
            path = self._project_root / self._template_specification_name(kind)
            spec = load_specification(path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Cannot load template project", str(error))
            return
        self._spec = spec
        self.spec_path_edit.setText(str(path))
        self.output_path_edit.setText(str(self._default_output_path(kind)))
        self._rebuild_fields()
        self._render_preview()

    def _choose_specification(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open project specification",
            str(self._workspace_paths().mesh_layout_projects),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        source_path = Path(path).expanduser().resolve()
        try:
            spec = load_specification(source_path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            QMessageBox.warning(self, "Invalid specification", str(error))
            return
        kind = spec["type"]
        template_kind_by_type = {
            "rectangle": "rectangle",
            "rectangle_vertical_midpoints": "rectangle_vertical_midpoints",
            "flat_plate": "flat_plate",
            "flat_plate_layered": "laminar_flat_plate_10_edges",
        }
        template_kind = template_kind_by_type.get(kind, kind)
        if source_path.parent.name == "laminar_flat_plate":
            template_kind = "laminar_flat_plate"
        elif source_path.parent.name == "laminar_flat_plate_10_edges":
            template_kind = "laminar_flat_plate_10_edges"
        index = self.template_combo.findData(template_kind)
        if index < 0:
            QMessageBox.warning(self, "Unsupported project", f"Unsupported type: {kind}")
            return
        try:
            source_path.relative_to(self._repository_root)
            is_in_repository = True
        except ValueError:
            is_in_repository = False
        if is_in_repository:
            canonical_template = self._template_project_source(template_kind) / self._template_specification_name(template_kind)
            if source_path != canonical_template.resolve():
                QMessageBox.warning(
                    self,
                    "Source tree is read-only",
                    "Copy this project outside the OneFLOW source tree before editing it.",
                )
                return
            project_root = ensure_template_project(
                self._repository_root,
                self._workspace_paths().mesh_layout_projects,
                template_kind,
            )
            source_path = project_root / self._template_specification_name(template_kind)
            spec = load_specification(source_path)
        else:
            project_root = source_path.parent
        self.template_combo.blockSignals(True)
        self.template_combo.setCurrentIndex(index)
        self.template_combo.blockSignals(False)
        self._project_root = project_root
        self._spec = spec
        self.spec_path_edit.setText(str(source_path))
        self.output_path_edit.setText(str(self._default_output_path(template_kind, project_root)))
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

        if self._spec["type"] in ("rectangle", "rectangle_vertical_midpoints"):
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
            if self._spec["type"] == "rectangle":
                for edge_index in range(4):
                    fields.append((("edgePointCounts", edge_index), f"Line {edge_index + 1} point count"))
                for line_id in range(1, 5):
                    fields.append((("boundaries", f"line{line_id}"), f"Boundary: Line {line_id}"))
            else:
                segment_names = ("bottom", "right lower", "right upper", "top", "left upper", "left lower")
                for segment_index, segment_name in enumerate(segment_names):
                    fields.append((("segmentPointCounts", segment_index), f"{segment_name.title()} segment point count"))
                for side in ("bottom", "right", "top", "left"):
                    fields.append((("boundaries", side), f"Boundary: {side.title()} side"))
            fields.append((("spacing",), "Spacing type"))
        else:
            fields = [
                (("upstreamDistance",), "Upstream distance"),
                (("plateLength",), "Plate length"),
                (("domainHeight",), "Domain height"),
                (("upstreamPointCount",), "Upstream line point count"),
                (("platePointCount",), "Plate line point count"),
                (("normalPointCount",), "Normal line point count"),
            ]
            distributions = self._spec.get("distributions", {})
            for distribution in ("upstream", "plate", "vertical", "topUpstream", "topPlate", "top"):
                if distribution not in distributions:
                    continue
                definition = distributions[distribution]
                if isinstance(definition, list):
                    for index, title in enumerate(("start spacing", "end spacing")):
                        fields.append((
                            ("distributions", distribution, index),
                            f"{distribution} {title}",
                        ))
                elif isinstance(definition, dict):
                    distribution_type = definition.get("type")
                    if distribution_type == "boundary_layer":
                        fields.extend([
                            (("distributions", distribution, "nearWallIntervalCount"),
                             f"{distribution} near-wall interval count"),
                            (("distributions", distribution, "nearWallDistance"),
                             f"{distribution} near-wall distance"),
                        ])
                    else:
                        for endpoint, title in (
                            ("startSpacing", "start spacing"),
                            ("endSpacing", "end spacing"),
                        ):
                            fields.append((
                                ("distributions", distribution, endpoint),
                                f"{distribution} {title}",
                            ))
            for boundary in (
                "upstreamBottom", "plate", "left", "interface", "right",
                "topUpstream", "topPlate",
            ):
                fields.append((("boundaries", boundary), f"Boundary: {boundary}"))

        for path, label in fields:
            value = self._get_path(self._spec, path)
            if path[0] == "boundaries":
                widget = QComboBox()
                for name, code, description in BOUNDARY_OPTIONS:
                    widget.addItem(name.title(), name)
                    widget.setItemData(
                        widget.count() - 1, description, Qt.ItemDataRole.ToolTipRole
                    )
                    widget.setItemData(
                        widget.count() - 1, code, Qt.ItemDataRole.UserRole + 1
                    )
                if isinstance(value, int) and not isinstance(value, bool):
                    selected_name = BOUNDARY_NAMES_BY_CODE.get(value)
                    if selected_name is None:
                        widget.addItem(
                            f"Custom code ({value}) - non-standard boundary code",
                            value,
                        )
                        selected_value = value
                    else:
                        selected_value = selected_name
                else:
                    selected_value = str(value).lower()
                selected_index = widget.findData(selected_value)
                if selected_index < 0:
                    widget.addItem(
                        f"Unknown value ({value}) - choose a supported boundary",
                        str(value),
                    )
                    selected_index = widget.count() - 1
                widget.setCurrentIndex(selected_index)
                widget.setToolTip(
                    "Choose a boundary type. Hover over an option for its explanation."
                )
            else:
                widget = QLineEdit(str(value))
                widget.setToolTip("Use 'auto' for an automatically selected endpoint spacing.")
            self._field_widgets[path] = widget
            self.fields_form.addRow(label, widget)

        self.reset_boundary_button.setEnabled(self._spec["type"] in ("flat_plate", "flat_plate_layered"))
        self._refresh_boundary_labels(self.show_boundary_codes.isChecked())

    def _refresh_boundary_labels(self, show_codes: bool) -> None:
        for path, widget in self._field_widgets.items():
            if path[0] != "boundaries" or not isinstance(widget, QComboBox):
                continue
            for index, (name, code, _) in enumerate(BOUNDARY_OPTIONS):
                if widget.itemData(index) == name:
                    label = name.title()
                    if show_codes:
                        label = f"{label} ({code})"
                    widget.setItemText(index, label)
            current_data = widget.currentData()
            if isinstance(current_data, int) and current_data not in BOUNDARY_NAMES_BY_CODE:
                widget.setItemText(
                    widget.findData(current_data),
                    f"Custom code ({current_data})" if show_codes else str(current_data),
                )

    def _reset_flat_plate_boundary_defaults(self) -> None:
        defaults = {
            "upstreamBottom": "symmetry",
            "plate": "wall",
            "left": "farfield",
            "interface": "interface",
            "right": "outflow",
            "topUpstream": "farfield",
            "topPlate": "farfield",
        }
        for name, value in defaults.items():
            widget = self._field_widgets.get(("boundaries", name))
            if isinstance(widget, QComboBox):
                index = widget.findData(value)
                if index >= 0:
                    widget.setCurrentIndex(index)

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
            original = self._get_path(spec, path)
            if isinstance(widget, QComboBox):
                value = widget.currentData()
                self._set_path(spec, path, value)
                continue
            raw = widget.text().strip()
            if path[0] == "corners":
                value = float(raw)
            elif path[0] in ("edgePointCounts", "segmentPointCounts") or (isinstance(path[-1], str) and path[-1].endswith("PointCount")):
                value = int(raw)
            elif path[0] == "distributions":
                if path[-1] == "nearWallIntervalCount":
                    value = int(raw)
                else:
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

    def _validated_fields(self) -> dict[str, Any] | None:
        try:
            spec = self._read_fields()
            compile_specification(spec)
            return spec
        except (ValueError, TypeError, KeyError) as error:
            QMessageBox.warning(self, "Invalid layout parameters", str(error))
            return None

    def _save_specification(self) -> None:
        spec = self._validated_fields()
        if spec is None:
            return
        current_text = self.spec_path_edit.text().strip()
        if not current_text:
            QMessageBox.warning(self, "Specification path required", "Choose a project JSON path.")
            return
        try:
            output_path = self._require_outside_repository(Path(current_text))
            save_specification(output_path, spec)
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.critical(self, "Cannot save specification", str(error))
            return
        self._spec = spec
        self._project_root = output_path.parent
        self.spec_path_edit.setText(str(output_path))
        self._rebuild_fields()
        self._render_preview()
        QMessageBox.information(self, "Specification saved", f"Saved project JSON to:\n{output_path}")

    def _save_specification_as(self) -> None:
        spec = self._validated_fields()
        if spec is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Save specification as",
            str(self._default_specification_path(self.template_combo.currentData())),
            "JSON files (*.json);;All files (*)",
        )
        if not path:
            return
        try:
            output_path = self._require_outside_repository(Path(path))
            save_specification(output_path, spec)
        except (OSError, ValueError, TypeError) as error:
            QMessageBox.critical(self, "Cannot save specification", str(error))
            return
        self._spec = spec
        self._project_root = output_path.parent
        self.spec_path_edit.setText(str(output_path))
        self.output_path_edit.setText(
            str(self._default_output_path(self.template_combo.currentData(), self._project_root))
        )
        self._rebuild_fields()
        self._render_preview()

    def _compile_layout(self) -> None:
        try:
            spec = self._read_fields()
            layout, svg = compile_specification(spec)
            if not str(self.output_path_edit.text()).strip():
                raise ValueError("Choose a layout output path.")
            output_path = self._require_outside_repository(
                Path(self.output_path_edit.text().strip())
            )
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
