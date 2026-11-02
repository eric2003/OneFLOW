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
    QFrame,
    QScrollArea,
    QSlider,
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
        self._latex_pixmap = QPixmap()
        self.latex_preview = QLabel("Render preview (mathtext)")
        self.latex_preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.latex_preview.setStyleSheet("background:#f7f7f7;")
        self.latex_scroll = QScrollArea()
        self.latex_scroll.setWidgetResizable(False)
        self.latex_scroll.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.latex_scroll.setMinimumHeight(180)
        self.latex_scroll.setMaximumHeight(240)
        self.latex_scroll.setFrameShape(QFrame.Shape.StyledPanel)
        self.latex_scroll.setStyleSheet("QScrollArea { background:#f7f7f7; }")
        self.latex_scroll.setWidget(self.latex_preview)
        self.latex_preview.setCursor(Qt.CursorShape.PointingHandCursor)
        self.latex_preview.mousePressEvent = self._on_latex_preview_click  # type: ignore[method-assign]

        latex_page = QWidget(); ll = QVBoxLayout(latex_page)
        ll.addWidget(QLabel("Plain text")); ll.addWidget(self.formula_plain, 1)
        ll.addWidget(QLabel("LaTeX source")); ll.addWidget(self.formula_latex, 1)
        lb = QHBoxLayout()
        b1 = QPushButton("Render preview"); b1.clicked.connect(self._render_latex)
        b2 = QPushButton("Copy LaTeX"); b2.clicked.connect(
            lambda: QApplication.clipboard().setText(self.formula_latex.toPlainText())
        )
        b_popup = QPushButton("Open large view...")
        b_popup.clicked.connect(self._open_latex_popup)
        lb.addWidget(b1); lb.addWidget(b2); lb.addWidget(b_popup); lb.addStretch(1)
        ll.addLayout(lb)
        ll.addWidget(QLabel("Preview (fixed frame; Open large view for zoom/pan)"))
        ll.addWidget(self.latex_scroll)

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
        if not e:
            return
        try:
            from verification.latex_preview import render_latex_preview_png
            preview_dir = UserDataPaths.discover().contract_workbench_preview
            preview_dir.mkdir(parents=True, exist_ok=True)
            out = preview_dir / f"{e.id}_latex.png"
            render_latex_preview_png(_latex(e), out, dpi=220, fontsize=20)
            self._latex_pixmap = QPixmap(str(out))
            vp = self.latex_scroll.viewport().size()
            if vp.width() < 40:
                vp.setWidth(520)
            if vp.height() < 40:
                vp.setHeight(180)
            fitted = self._latex_pixmap.scaled(
                vp,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            self.latex_preview.setPixmap(fitted)
            self.latex_preview.resize(fitted.size())
            self.latex_preview.setToolTip(
                "Click to open large view (zoom toward center / pan / reset)."
            )
        except Exception as exc:
            self._latex_pixmap = QPixmap()
            self.latex_preview.setText(
                f"Preview unavailable ({exc}). Use LaTeX source above."
            )
            self.latex_preview.adjustSize()

    def _on_latex_preview_click(self, event) -> None:  # noqa: ANN001
        if event.button() == Qt.MouseButton.LeftButton:
            self._open_latex_popup()
            return
        QLabel.mousePressEvent(self.latex_preview, event)

    def _open_latex_popup(self) -> None:
        from PySide6.QtCore import QPoint, QTimer

        if self._latex_pixmap is None or self._latex_pixmap.isNull():
            self._render_latex()
        if self._latex_pixmap is None or self._latex_pixmap.isNull():
            QMessageBox.information(
                self, "LaTeX preview", "Render a preview first."
            )
            return

        dlg = QDialog(self)
        dlg.setWindowTitle("Formula large view")
        dlg.setWindowFlags(
            Qt.WindowType.Window
            | Qt.WindowType.WindowTitleHint
            | Qt.WindowType.WindowSystemMenuHint
            | Qt.WindowType.WindowMinMaxButtonsHint
            | Qt.WindowType.WindowCloseButtonHint
        )
        dlg.resize(920, 560)

        lay = QVBoxLayout(dlg)
        scroll = QScrollArea()
        scroll.setWidgetResizable(False)
        scroll.setAlignment(Qt.AlignmentFlag.AlignCenter)
        scroll.setStyleSheet("QScrollArea { background:#ffffff; }")
        img = QLabel()
        img.setAlignment(Qt.AlignmentFlag.AlignCenter)
        img.setStyleSheet("background:#ffffff;")
        scroll.setWidget(img)
        lay.addWidget(scroll, 1)

        state = {"z": 1.0}

        def _viewport_center_on_content() -> tuple[float, float]:
            """Content coordinates currently under the viewport center."""
            vp = scroll.viewport().size()
            hbar = scroll.horizontalScrollBar()
            vbar = scroll.verticalScrollBar()
            cx = hbar.value() + vp.width() / 2.0
            cy = vbar.value() + vp.height() / 2.0
            return cx, cy

        def _scroll_so_content_point_at_center(cx: float, cy: float) -> None:
            """Keep (cx, cy) in content coords under the viewport center after zoom."""
            vp = scroll.viewport().size()
            hbar = scroll.horizontalScrollBar()
            vbar = scroll.verticalScrollBar()
            hbar.setValue(int(round(cx - vp.width() / 2.0)))
            vbar.setValue(int(round(cy - vp.height() / 2.0)))

        def _show(*, preserve_center: bool = False) -> None:
            old_cx = old_cy = 0.0
            old_w = max(1, img.width())
            old_h = max(1, img.height())
            if preserve_center and not img.pixmap().isNull():
                old_cx, old_cy = _viewport_center_on_content()

            z = state["z"]
            w = max(1, int(self._latex_pixmap.width() * z))
            h = max(1, int(self._latex_pixmap.height() * z))
            pix = self._latex_pixmap.scaled(
                w,
                h,
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
            img.setPixmap(pix)
            img.resize(pix.size())
            zlab.setText(f"Zoom {int(round(z * 100))}%")

            if preserve_center:
                # Map previous center through scale change.
                scale_x = pix.width() / float(old_w)
                scale_y = pix.height() / float(old_h)
                _scroll_so_content_point_at_center(old_cx * scale_x, old_cy * scale_y)
            else:
                # Fit / reset: center content in the viewport when possible.
                vp = scroll.viewport().size()
                if pix.width() <= vp.width():
                    scroll.horizontalScrollBar().setValue(0)
                else:
                    scroll.horizontalScrollBar().setValue(
                        max(0, (pix.width() - vp.width()) // 2)
                    )
                if pix.height() <= vp.height():
                    scroll.verticalScrollBar().setValue(0)
                else:
                    scroll.verticalScrollBar().setValue(
                        max(0, (pix.height() - vp.height()) // 2)
                    )

        def zin() -> None:
            state["z"] = min(5.0, state["z"] * 1.25)
            _show(preserve_center=True)

        def zout() -> None:
            state["z"] = max(0.25, state["z"] / 1.25)
            _show(preserve_center=True)

        def zreset() -> None:
            vp_w = max(100, scroll.viewport().width() - 16)
            if self._latex_pixmap.width() > 0:
                state["z"] = min(1.0, vp_w / float(self._latex_pixmap.width()))
            else:
                state["z"] = 1.0
            _show(preserve_center=False)

        def pan(dx: int, dy: int) -> None:
            scroll.horizontalScrollBar().setValue(
                scroll.horizontalScrollBar().value() + dx
            )
            scroll.verticalScrollBar().setValue(
                scroll.verticalScrollBar().value() + dy
            )

        bar = QHBoxLayout()
        b_out = QPushButton("Zoom -")
        b_out.clicked.connect(zout)
        b_in = QPushButton("Zoom +")
        b_in.clicked.connect(zin)
        b_rs = QPushButton("Reset / fit width")
        b_rs.clicked.connect(zreset)
        b_l = QPushButton("Left")
        b_l.clicked.connect(lambda: pan(-80, 0))
        b_r = QPushButton("Right")
        b_r.clicked.connect(lambda: pan(80, 0))
        b_u = QPushButton("Up")
        b_u.clicked.connect(lambda: pan(0, -80))
        b_d = QPushButton("Down")
        b_d.clicked.connect(lambda: pan(0, 80))
        zlab = QLabel()
        for w in (b_out, b_in, b_rs, b_l, b_r, b_u, b_d, zlab):
            bar.addWidget(w)
        bar.addStretch(1)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.accept)
        bar.addWidget(close_btn)
        lay.addLayout(bar)
        hint = QLabel(
            "Window opens centered. Zoom keeps the viewport center fixed "
            "(formula stays centered while scaling). Pan moves the view; "
            "Reset / fit width shows the full line width."
        )
        hint.setWordWrap(True)
        lay.addWidget(hint)

        def _center_and_fit() -> None:
            parent = self.window()
            if parent is not None:
                top_left = parent.mapToGlobal(QPoint(0, 0))
                x = top_left.x() + (parent.width() - dlg.width()) // 2
                y = top_left.y() + (parent.height() - dlg.height()) // 2
                dlg.move(max(0, x), max(0, y))
            else:
                screen = dlg.screen()
                if screen is not None:
                    ag = screen.availableGeometry()
                    dlg.move(
                        ag.center().x() - dlg.width() // 2,
                        ag.center().y() - dlg.height() // 2,
                    )
            zreset()

        dlg.show()
        QTimer.singleShot(0, _center_and_fit)
        dlg.exec()

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
