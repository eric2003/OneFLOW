"""Contract Workbench: browse formula knowledge + theory curves.

Developer aid: inspect formula, principle, background, references, pure APIs,
and optional static plots. Paths are resolved from repository_root only.
"""

from __future__ import annotations

import sys
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSplitter,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)


def _ensure_tools_on_path(repository_root: Path) -> None:
    tools_root = (repository_root / "tools").resolve()
    tools_str = str(tools_root)
    if tools_str not in sys.path:
        sys.path.insert(0, tools_str)


class VerificationCurveDialog(QDialog):
    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        _ensure_tools_on_path(self._repository_root)

        from verification.curve_catalog import list_theory_entries

        self._entries = list_theory_entries()
        self._id_to_entry = {e.id: e for e in self._entries}
        self._preview_png: Path | None = None

        self.setWindowTitle("Contract Workbench — Formula & Theory Curves")
        self.resize(1100, 720)

        root = QVBoxLayout(self)
        intro = QLabel(
            "Browse numerical knowledge by taxonomy (family / model / section). "
            "Each item shows formula, principle, background, references, and "
            "linked pure APIs. Plots use repository-relative theory CSV files "
            "(portable across machines). Not a substitute for GoogleTest or "
            "case regression."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)

        splitter = QSplitter()
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Taxonomy / entry"])
        self.tree.setMinimumWidth(280)
        self._fill_tree()
        self.tree.currentItemChanged.connect(self._on_tree)
        splitter.addWidget(self.tree)

        right = QWidget()
        right_l = QVBoxLayout(right)
        self.tabs = QTabWidget()
        self.overview = QPlainTextEdit()
        self.overview.setReadOnly(True)
        self.formula = QPlainTextEdit()
        self.formula.setReadOnly(True)
        self.formula.setStyleSheet("font-family: Consolas, monospace;")
        self.principle = QPlainTextEdit()
        self.principle.setReadOnly(True)
        self.refs = QPlainTextEdit()
        self.refs.setReadOnly(True)
        self.plot_panel = QWidget()
        plot_l = QVBoxLayout(self.plot_panel)
        self.plot_label = QLabel("Select an entry with a theory curve to plot.")
        self.plot_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.plot_label.setMinimumHeight(320)
        self.plot_label.setStyleSheet("background: #1e1e1e; color: #ddd;")
        plot_l.addWidget(self.plot_label, 1)
        plot_btns = QHBoxLayout()
        self.plot_button = QPushButton("Render plot")
        self.plot_button.clicked.connect(self._render_plot)
        self.export_button = QPushButton("Export PNG...")
        self.export_button.clicked.connect(self._export_png)
        plot_btns.addWidget(self.plot_button)
        plot_btns.addWidget(self.export_button)
        plot_btns.addStretch(1)
        plot_l.addLayout(plot_btns)

        self.tabs.addTab(self.overview, "Overview")
        self.tabs.addTab(self.formula, "Formula")
        self.tabs.addTab(self.principle, "Principle & background")
        self.tabs.addTab(self.refs, "References")
        self.tabs.addTab(self.plot_panel, "Plot")
        right_l.addWidget(self.tabs, 1)

        close_row = QHBoxLayout()
        close_row.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.accept)
        close_row.addWidget(close_btn)
        right_l.addLayout(close_row)

        splitter.addWidget(right)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter, 1)

        # Select first leaf if any
        self._select_first_leaf()

    def _fill_tree(self) -> None:
        self.tree.clear()
        # family -> genus -> section -> entries
        nodes: dict[tuple, QTreeWidgetItem] = {}
        for e in self._entries:
            fam_key = (e.family,)
            if fam_key not in nodes:
                nodes[fam_key] = QTreeWidgetItem([e.family])
                self.tree.addTopLevelItem(nodes[fam_key])
            gen_key = (e.family, e.genus)
            if gen_key not in nodes:
                nodes[gen_key] = QTreeWidgetItem([e.genus])
                nodes[fam_key].addChild(nodes[gen_key])
            sec_key = (e.family, e.genus, e.section)
            if sec_key not in nodes:
                nodes[sec_key] = QTreeWidgetItem([e.section])
                nodes[gen_key].addChild(nodes[sec_key])
            leaf = QTreeWidgetItem([e.title])
            leaf.setData(0, Qt.ItemDataRole.UserRole, e.id)
            nodes[sec_key].addChild(leaf)
        self.tree.expandAll()

    def _select_first_leaf(self) -> None:
        def walk(item: QTreeWidgetItem) -> QTreeWidgetItem | None:
            if item.data(0, Qt.ItemDataRole.UserRole):
                return item
            for i in range(item.childCount()):
                found = walk(item.child(i))
                if found:
                    return found
            return None

        for i in range(self.tree.topLevelItemCount()):
            found = walk(self.tree.topLevelItem(i))
            if found:
                self.tree.setCurrentItem(found)
                return

    def _current_entry(self):
        item = self.tree.currentItem()
        if item is None:
            return None
        entry_id = item.data(0, Qt.ItemDataRole.UserRole)
        if not entry_id:
            return None
        return self._id_to_entry.get(entry_id)

    def _on_tree(self, _current, _previous) -> None:
        e = self._current_entry()
        self._preview_png = None
        self.plot_label.setPixmap(QPixmap())
        if e is None:
            self.overview.setPlainText("Select a leaf entry in the tree.")
            self.formula.clear()
            self.principle.clear()
            self.refs.clear()
            self.plot_label.setText("Select an entry with a theory curve to plot.")
            self.plot_button.setEnabled(False)
            self.export_button.setEnabled(False)
            return

        from verification.curve_catalog import csv_path, doc_path

        csv = csv_path(self._repository_root, e)
        doc = doc_path(self._repository_root, e)
        tags = ", ".join(e.role_tags)
        apis = ", ".join(e.pure_apis)
        self.overview.setPlainText(
            f"id: {e.id}\n"
            f"taxonomy: {e.family} / {e.genus} / {e.section}\n"
            f"title: {e.title}\n"
            f"role tags: {tags}\n\n"
            f"pure APIs: {apis}\n"
            f"doc (repo-relative): {e.doc_relpath or '(none)'}\n"
            f"doc exists: {bool(doc and doc.is_file())}\n"
            f"theory CSV (repo-relative): {e.csv_relpath or '(none)'}\n"
            f"CSV exists: {bool(csv and csv.is_file())}\n\n"
            f"Proves: {e.proves}\n"
            f"Does not prove: {e.does_not_prove}\n"
        )
        self.formula.setPlainText(e.formula)
        self.principle.setPlainText(
            f"Principle\n---------\n{e.principle}\n\n"
            f"Background\n----------\n{e.background}\n"
        )
        self.refs.setPlainText("\n".join(f"- {r}" for r in e.references))
        has_plot = bool(e.plot_module and e.csv_relpath)
        self.plot_button.setEnabled(has_plot)
        self.export_button.setEnabled(has_plot)
        if has_plot:
            self.plot_label.setText(
                "Click Render plot to draw from the theory CSV "
                "(path resolved under this repository root)."
            )
        else:
            self.plot_label.setText(
                "This entry is formula knowledge only (no theory curve file)."
            )

    def _render_plot(self) -> None:
        e = self._current_entry()
        if e is None:
            return
        try:
            from verification.curve_catalog import render_png

            # Portable temp under repo tools/verification/data/_preview
            preview_dir = (
                self._repository_root / "tools" / "verification" / "data" / "_preview"
            )
            preview_dir.mkdir(parents=True, exist_ok=True)
            out = preview_dir / f"{e.id}_preview.png"
            render_png(self._repository_root, e, out)
            self._preview_png = out
            pix = QPixmap(str(out))
            if pix.isNull():
                raise RuntimeError("failed to load rendered PNG")
            self.plot_label.setPixmap(
                pix.scaled(
                    self.plot_label.size(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
            )
        except Exception as exc:
            QMessageBox.warning(self, "Plot failed", str(exc))

    def _export_png(self) -> None:
        e = self._current_entry()
        if e is None:
            return
        out, _ = QFileDialog.getSaveFileName(
            self, "Export theory curve PNG", f"{e.id}.png", "PNG (*.png)"
        )
        if not out:
            return
        try:
            from verification.curve_catalog import render_png

            render_png(self._repository_root, e, Path(out))
            QMessageBox.information(self, "Export", f"Wrote\n{out}")
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))
