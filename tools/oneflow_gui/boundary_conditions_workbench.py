"""Dedicated Boundary Conditions Workbench."""
from __future__ import annotations
import json
from pathlib import Path
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QVBoxLayout, QLabel, QLineEdit, QListWidget, QListWidgetItem, QSplitter, QTabWidget, QPlainTextEdit, QTableWidget, QTableWidgetItem, QWidget

REGISTRY = Path("tools/verification/boundary_conditions/registry.json")

def load_boundary_registry(root: Path) -> dict:
    data = json.loads((Path(root) / REGISTRY).read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError("Unsupported boundary registry schema")
    entries = data.get("conditions", [])
    seen = set()
    for item in entries:
        required = ("id", "name", "formulation", "parameters", "mappings", "implementation", "verification")
        if not all(key in item for key in required):
            raise ValueError("Incomplete boundary registry entry")
        if item["id"] in seen:
            raise ValueError("Duplicate boundary condition id: " + item["id"])
        seen.add(item["id"])
    if not entries:
        raise ValueError("Boundary registry is empty")
    return data

class BoundaryConditionsWorkbenchDialog(QDialog):
    def __init__(self, repository_root: Path, parent=None):
        super().__init__(parent)
        self.registry = load_boundary_registry(repository_root)
        self.entries = tuple(self.registry["conditions"])
        self.by_id = {entry["id"]: entry for entry in self.entries}
        self.setWindowTitle("Boundary Conditions Workbench")
        self.resize(1280, 800)
        root = QVBoxLayout(self)
        root.addWidget(QLabel("Boundary Conditions Workbench"))
        intro = QLabel("Physical formulations, equations, parameters, format mappings, OneFLOW implementation status and verification evidence. Mapping does not imply solver support.")
        intro.setWordWrap(True)
        root.addWidget(intro)
        self.search_edit = QLineEdit()
        self.search_edit.setObjectName("boundarySearch")
        self.search_edit.setPlaceholderText("Search semantic IDs, CFL3D types, parameters or mappings...")
        self.search_edit.textChanged.connect(self._populate)
        root.addWidget(self.search_edit)
        split = QSplitter()
        self.condition_list = QListWidget()
        self.condition_list.currentItemChanged.connect(self._select)
        split.addWidget(self.condition_list)
        detail = QWidget()
        from PySide6.QtWidgets import QVBoxLayout as VBox
        layout = VBox(detail)
        self.heading = QLabel("Select a boundary condition")
        layout.addWidget(self.heading)
        self.tabs = QTabWidget()
        self.views = {}
        for name in ("Overview", "Formulation", "Parameters", "Mappings", "Implementation", "Verification", "References"):
            if name == "Mappings":
                view = QTableWidget(0, 3)
                view.setHorizontalHeaderLabels(["System", "Type / Name", "Relation"])
                view.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
                view.horizontalHeader().setStretchLastSection(True)
            else:
                view = QPlainTextEdit()
                view.setReadOnly(True)
            self.views[name] = view
            self.tabs.addTab(view, name)
        layout.addWidget(self.tabs)
        split.addWidget(detail)
        split.setStretchFactor(1, 3)
        root.addWidget(split, 1)
        self._populate()
        if self.condition_list.count():
            self.condition_list.setCurrentRow(0)

    def _populate(self, query=""):
        old = self.condition_list.currentItem()
        old_id = old.data(Qt.ItemDataRole.UserRole) if old else None
        self.condition_list.clear()
        for entry in self.entries:
            if query and query.casefold() not in json.dumps(entry).casefold():
                continue
            row = QListWidgetItem(entry["name"])
            row.setData(Qt.ItemDataRole.UserRole, entry["id"])
            self.condition_list.addItem(row)
        selected = next((i for i in range(self.condition_list.count()) if self.condition_list.item(i).data(Qt.ItemDataRole.UserRole) == old_id), 0)
        if self.condition_list.count():
            self.condition_list.setCurrentRow(selected)

    def _text(self, tab, lines):
        self.views[tab].setPlainText("\n".join(lines))

    def _select(self, current, _previous=None):
        if current is None:
            return
        entry = self.by_id[current.data(Qt.ItemDataRole.UserRole)]
        cfl = entry["cfl3d"]
        self.heading.setText(entry["name"])
        self._text("Overview", [entry["description"], "", "CFL3D " + str(cfl["version"]) + " type " + str(cfl["bctype"]), "Manual section " + cfl["section"], "", "Data location: " + entry["data_location"]])
        self._text("Formulation", ["Equations / relations:", *entry["formulation"]["equations"], "", entry["formulation"].get("notes", "")])
        self._text("Parameters", [item["name"] + " — required" for item in entry["parameters"]] or ["No additional parameters declared."])
        mapping = entry["mappings"]
        rows = [["CFL3D", str(cfl["bctype"]), cfl["section"]], ["CGNS", str(mapping["cgns"]), "candidate mapping"], ["OneFLOW", str(mapping["oneflow"]), entry["implementation"]["status"]]]
        table = self.views["Mappings"]
        table.setRowCount(len(rows))
        for r, values in enumerate(rows):
            for c, value in enumerate(values):
                table.setItem(r, c, QTableWidgetItem(str(value)))
        self._text("Implementation", ["Status: " + entry["implementation"]["status"], entry["implementation"]["numerical_realization"]])
        self._text("Verification", [str(k) + ": " + str(v) for k, v in entry["verification"].items()])
        source = self.registry["source"]
        self._text("References", [source["name"], source["url"], source["document"], "Chapter " + source["chapter"] + ", Section " + cfl["section"]])