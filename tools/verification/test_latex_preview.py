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

    def test_parameter_contract_latex_uses_single_command_escapes(self) -> None:
        import json

        contracts_root = (
            Path(__file__).resolve().parent
            / "contracts"
            / "parameter"
        )
        for contract_name in ("inflow_intensity", "inflow_viscosity", "rprod"):
            with self.subTest(contract=contract_name):
                definition = json.loads(
                    (contracts_root / f"{contract_name}.json").read_text(
                        encoding="utf-8"
                    )
                )
                formulas = definition["knowledge"]["formulas"]
                self.assertTrue(formulas)
                if contract_name == "inflow_intensity":
                    self.assertEqual(
                        formulas[0]["latex"],
                        r"k_e = \frac{3}{2} (0.01 I)^2",
                    )
                for formula in formulas:
                    latex = formula["latex"]
                    self.assertNotIn(
                        "\\\\", latex,
                        "LaTeX commands must start with one backslash.",
                    )
                    with tempfile.TemporaryDirectory() as td:
                        output = Path(td) / "formula.png"
                        render_latex_preview_png(latex, output)
                        self.assertTrue(output.is_file())
                        self.assertGreater(output.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
