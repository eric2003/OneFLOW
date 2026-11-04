from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from .latex_preview import render_latex_preview_png, sanitize_for_mathtext


class LatexPreviewTest(unittest.TestCase):
    def test_sanitize_strips_big(self) -> None:
        s = r"\mathrm{arg}_1=\min\Big(\max(\cdots),\cdots\Big)"
        out = sanitize_for_mathtext(s)
        self.assertNotIn(r"\Big", out)
        self.assertNotIn(r"\mathrm", out)

    def test_sanitize_preserves_latex_spacing_commands(self) -> None:
        source = r"x \quad y + \frac{a}{b}"
        sanitized = sanitize_for_mathtext(source)
        self.assertIn(r"\quad", sanitized)
        self.assertIn(r"\frac", sanitized)

    def test_render_does_not_raise_on_full_latex(self) -> None:
        body = r"""
F_1=\tanh(arg_1^4)
\mathsf{arg}_1=\min\Big(\max(\cdots),\,\cdots\Big)
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
        for contract_path in sorted(contracts_root.glob("*.json")):
            contract_name = contract_path.stem
            with self.subTest(contract=contract_name):
                definition = json.loads(
                    contract_path.read_text(encoding="utf-8")
                )
                formulas = definition.get("knowledge", {}).get("formulas", [])
                for formula in formulas:
                    latex = formula["latex"]
                    if contract_name == "inflow_intensity":
                        self.assertEqual(
                            latex,
                            r"k_e = \frac{3}{2} (0.01 I)^2",
                        )
                    self.assertNotIn(
                        chr(92) * 2,
                        latex,
                        "LaTeX commands must start with one backslash.",
                    )
                    with tempfile.TemporaryDirectory() as td:
                        output = Path(td) / "formula.png"
                        render_latex_preview_png(latex, output)
                        self.assertTrue(output.is_file())
                        self.assertGreater(output.stat().st_size, 0)



    def test_time_step_strategy_formula_uses_balanced_multiline_labels(self) -> None:
        import json

        contract_path = (
            Path(__file__).resolve().parent
            / "contracts"
            / "parameter"
            / "time_step_strategy.json"
        )
        definition = json.loads(contract_path.read_text(encoding="utf-8"))
        formula = definition["knowledge"]["formulas"][0]
        latex = formula["latex"]
        self.assertEqual(
            latex.splitlines(),
            [
                r"S(m)=\mathsf{local},\quad m=0",
                r"S(m)=\mathsf{global},\quad m=1",
                r"S(m)=\mathsf{local-global},\quad \mathsf{otherwise}",
            ],
        )
        self.assertNotIn(r"\text{-}", latex)
        with tempfile.TemporaryDirectory() as td:
            output = Path(td) / "time-step-strategy.png"
            render_latex_preview_png(latex, output, fontsize=26)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
