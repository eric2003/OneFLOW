"""Contract Workbench knowledge browser."""

from __future__ import annotations

import sys
from pathlib import Path

from .user_data_paths import UserDataPaths

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
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
    s = str(tools_root)
    if s not in sys.path:
        sys.path.insert(0, s)


def _plain(e) -> str:
    return getattr(e, "formula_plain", None) or getattr(e, "formula", "") or ""


def _latex(e) -> str:
    return getattr(e, "formula_latex", None) or _plain(e)


def _sketch(e) -> str:
    return getattr(e, "code_sketch", None) or "// No gated sketch yet.\n"


def _sources(e) -> tuple:
    return getattr(e, "source_relpaths", ()) or ()


class VerificationCurveDialog(QDialog):
    def __init__(self, repository_root: Path, parent=None) -> None:
        super().__init__(parent)
        self._repository_root = repository_root.resolve()
        _ensure_tools_on_path(self._repository_root)
        from verification.curve_catalog import list_theory_entries
        self._entries = list_theory_entries()
        self._map = {e.id: e for e in self._entries}
        self.setWindowTitle("Contract Workbench — Formula knowledge")
        self.resize(1180, 760)
        root = QVBoxLayout(self)
        intro = QLabel(
            "Taxonomy: turbulence (SA/SST) and mean-flow (NS-LUSGS / NS-inviscid). "
            "LaTeX preview is mathtext best-effort; full LaTeX source is always available. "
            "Paths are repository-relative."
        )
        intro.setWordWrap(True)
        root.addWidget(intro)
        split = QSplitter()
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Taxonomy / entry"])
        self.tree.setMinimumWidth(300)
        self._fill_tree()
        self.tree.currentItemChanged.connect(self._on_tree)
        split.addWidget(self.tree)
        right = QWidget()
        rl = QVBoxLayout(right)
        self.tabs = QTabWidget()
        self.overview = QPlainTextEdit(); self.overview.setReadOnly(True)
        self.formula_plain = QPlainTextEdit(); self.formula_plain.setReadOnly(True)
        self.formula_plain.setStyleSheet("font-family: Consolas, monospace;")
        self.formula_latex = QPlainTextEdit(); self.formula_latex.setReadOnly(True)
        self.formula_latex.setStyleSheet("font-family: Consolas, monospace;")
        self.latex_preview = QLabel("Render preview (mathtext)")
        self.latex_preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.latex_preview.setMinimumHeight(180)
        self.latex_preview.setStyleSheet("background:#f7f7f7;")
        latex_page = QWidget(); ll = QVBoxLayout(latex_page)
        ll.addWidget(QLabel("Plain text")); ll.addWidget(self.formula_plain, 1)
        ll.addWidget(QLabel("LaTeX source")); ll.addWidget(self.formula_latex, 1)
        lb = QHBoxLayout()
        b1 = QPushButton("Render preview"); b1.clicked.connect(self._render_latex)
        b2 = QPushButton("Copy LaTeX"); b2.clicked.connect(lambda: QApplication.clipboard().setText(self.formula_latex.toPlainText()))
        lb.addWidget(b1); lb.addWidget(b2); lb.addStretch(1)
        ll.addLayout(lb); ll.addWidget(self.latex_preview)
        self.principle = QPlainTextEdit(); self.principle.setReadOnly(True)
        self.refs = QPlainTextEdit(); self.refs.setReadOnly(True)
        self.doc_view = QPlainTextEdit(); self.doc_view.setReadOnly(True)
        self.sketch_view = QPlainTextEdit(); self.sketch_view.setReadOnly(True)
        self.sketch_view.setStyleSheet("font-family: Consolas, monospace;")
        self.allow_sketch = QCheckBox("Gated sketch: I will not auto-write production")
        sk = QWidget(); skl = QVBoxLayout(sk); skl.addWidget(self.allow_sketch); skl.addWidget(self.sketch_view, 1)
        self.copy_sk = QPushButton("Copy sketch"); self.copy_sk.setEnabled(False)
        self.allow_sketch.toggled.connect(self.copy_sk.setEnabled)
        self.copy_sk.clicked.connect(lambda: self.allow_sketch.isChecked() and QApplication.clipboard().setText(self.sketch_view.toPlainText()))
        skl.addWidget(self.copy_sk)
        self.plot_label = QLabel("Curve plot"); self.plot_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.plot_label.setMinimumHeight(280); self.plot_label.setStyleSheet("background:#1e1e1e;color:#ddd;")
        pp = QWidget(); pl = QVBoxLayout(pp); pl.addWidget(self.plot_label, 1)
        pb = QHBoxLayout(); self.plot_btn = QPushButton("Render curve"); self.plot_btn.clicked.connect(self._render_curve)
        self.export_btn = QPushButton("Export PNG..."); self.export_btn.clicked.connect(self._export_curve)
        pb.addWidget(self.plot_btn); pb.addWidget(self.export_btn); pb.addStretch(1); pl.addLayout(pb)
        self.tabs.addTab(self.overview, "Overview")
        self.tabs.addTab(latex_page, "Formula")
        self.tabs.addTab(self.principle, "Principle")
        self.tabs.addTab(self.refs, "References")
        self.tabs.addTab(self.doc_view, "Sphinx excerpt")
        self.tabs.addTab(sk, "Code sketch")
        self.tabs.addTab(pp, "Plot")
        rl.addWidget(self.tabs, 1)
        nav = QHBoxLayout()
        od = QPushButton("Open Sphinx .rst..."); od.clicked.connect(self._open_doc)
        osrc = QPushButton("Open source files..."); osrc.clicked.connect(self._open_src)
        nav.addWidget(od); nav.addWidget(osrc); nav.addStretch(1)
        cl = QPushButton("Close"); cl.clicked.connect(self.accept); nav.addWidget(cl)
        rl.addLayout(nav)
        split.addWidget(right); split.setStretchFactor(1, 1); root.addWidget(split, 1)
        self._select_first()

    def _fill_tree(self) -> None:
        self.tree.clear(); nodes = {}
        for e in self._entries:
            fk, gk, sk = (e.family,), (e.family, e.genus), (e.family, e.genus, e.section)
            if fk not in nodes:
                nodes[fk] = QTreeWidgetItem([e.family]); self.tree.addTopLevelItem(nodes[fk])
            if gk not in nodes:
                nodes[gk] = QTreeWidgetItem([e.genus]); nodes[fk].addChild(nodes[gk])
            if sk not in nodes:
                nodes[sk] = QTreeWidgetItem([e.section]); nodes[gk].addChild(nodes[sk])
            leaf = QTreeWidgetItem([e.title]); leaf.setData(0, Qt.ItemDataRole.UserRole, e.id)
            nodes[sk].addChild(leaf)
        self.tree.expandAll()

    def _select_first(self) -> None:
        def walk(it):
            if it.data(0, Qt.ItemDataRole.UserRole): return it
            for i in range(it.childCount()):
                f = walk(it.child(i))
                if f: return f
            return None
        for i in range(self.tree.topLevelItemCount()):
            f = walk(self.tree.topLevelItem(i))
            if f: self.tree.setCurrentItem(f); return

    def _cur(self):
        it = self.tree.currentItem()
        if not it: return None
        eid = it.data(0, Qt.ItemDataRole.UserRole)
        return self._map.get(eid) if eid else None

    def _on_tree(self, *a) -> None:
        from verification.curve_catalog import csv_path, doc_path, read_doc_excerpt
        e = self._cur(); self.allow_sketch.setChecked(False)
        self.latex_preview.setPixmap(QPixmap()); self.plot_label.setPixmap(QPixmap())
        if e is None:
            self.overview.setPlainText("Select a leaf."); return
        csv = csv_path(self._repository_root, e)
        self.overview.setPlainText(
            f"id: {e.id}\ntaxonomy: {e.family} / {e.genus} / {e.section}\n"
            f"title: {e.title}\nrole tags: {', '.join(e.role_tags)}\n"
            f"pure APIs: {', '.join(e.pure_apis)}\ndoc: {e.doc_relpath}\n"
            f"sources:\n  " + "\n  ".join(_sources(e)) + "\n"
            f"CSV: {e.csv_relpath} exists={bool(csv and csv.is_file())}\n\n"
            f"Proves: {e.proves}\nDoes not prove: {e.does_not_prove}\n"
        )
        self.formula_plain.setPlainText(_plain(e))
        self.formula_latex.setPlainText(_latex(e))
        self.principle.setPlainText(f"Principle\n---------\n{e.principle}\n\nBackground\n----------\n{e.background}\n")
        self.refs.setPlainText("\n".join(f"- {r}" for r in e.references))
        self.doc_view.setPlainText(read_doc_excerpt(self._repository_root, e))
        self.sketch_view.setPlainText(_sketch(e))
        has = bool(e.plot_module and e.csv_relpath)
        self.plot_btn.setEnabled(has); self.export_btn.setEnabled(has)
        self.plot_label.setText("Click Render curve." if has else "No theory curve.")

    def _render_latex(self) -> None:
        e = self._cur()
        if not e: return
        try:
            from verification.latex_preview import render_latex_preview_png
            preview_dir = UserDataPaths.discover().contract_workbench_preview
            preview_dir.mkdir(parents=True, exist_ok=True)
            out = preview_dir / f"{e.id}_latex.png"
            render_latex_preview_png(_latex(e), out)
            pix = QPixmap(str(out))
            self.latex_preview.setPixmap(pix.scaled(self.latex_preview.size(), Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
        except Exception as exc:
            self.latex_preview.setText(f"Preview unavailable ({exc}). Use LaTeX source above.")

    def _render_curve(self) -> None:
        e = self._cur()
        if not e: return
        try:
            from verification.curve_catalog import render_png
            preview_dir = UserDataPaths.discover().contract_workbench_preview
            preview_dir.mkdir(parents=True, exist_ok=True)
            out = preview_dir / f"{e.id}_curve.png"
            render_png(self._repository_root, e, out)
            pix = QPixmap(str(out))
            self.plot_label.setPixmap(pix.scaled(self.plot_label.size(), Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
        except Exception as exc:
            QMessageBox.warning(self, "Plot failed", str(exc))

    def _export_curve(self) -> None:
        e = self._cur()
        if not e: return
        path, _ = QFileDialog.getSaveFileName(self, "Export PNG", f"{e.id}.png", "PNG (*.png)")
        if not path: return
        try:
            from verification.curve_catalog import render_png
            render_png(self._repository_root, e, Path(path))
            QMessageBox.information(self, "Export", f"Wrote\n{path}")
        except Exception as exc:
            QMessageBox.warning(self, "Export failed", str(exc))

    def _open_doc(self) -> None:
        from verification.curve_catalog import doc_path
        e = self._cur()
        if not e: return
        p = doc_path(self._repository_root, e)
        if not p or not p.is_file():
            QMessageBox.warning(self, "Doc", "Not found"); return
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(p)))

    def _open_src(self) -> None:
        from verification.curve_catalog import resolve_repo_path
        e = self._cur()
        if not e: return
        n = 0
        for rel in _sources(e):
            p = resolve_repo_path(self._repository_root, rel)
            if p and p.is_file():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(p))); n += 1
        if not n:
            QMessageBox.warning(self, "Source", "No files found")
