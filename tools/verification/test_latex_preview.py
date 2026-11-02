from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from latex_preview import render_latex_preview_png, sanitize_for_mathtext


class LatexPreviewTest(unittest.TestCase):
    def test_sanitize_strips_big(self) -> None:
        s = r"\mathrm{arg}_1=\min\Big(\max(\cdots),\cdots\Big)"
        out = sanitize_for_mathtext(s)
        self.assertNotIn(r"\Big", out)
        self.assertNotIn(r"\mathrm", out)

    def test_render_does_not_raise_on_full_latex(self) -> None:
        body = r"""
F_1=\tanh(arg_1^4)
\mathrm{arg}_1=\min\Big(\max(\cdots),\,\cdots\Big)
"""
        with tempfile.TemporaryDirectory() as td:
            out = Path(td) / "x.png"
            render_latex_preview_png(body, out)
            self.assertTrue(out.is_file())


if __name__ == "__main__":
    unittest.main()
