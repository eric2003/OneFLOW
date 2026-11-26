from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.verification.change_gate import select_test_commands


class ChangeGateSelectionTest(unittest.TestCase):
    def setUp(self):
        self.registry = (
            ("mesh geometry", ("-m", "unittest", "grid.tests.python.stl_to_gmsh.test_geometry")),
            ("surface import", ("-m", "unittest", "grid.tests.python.stl_import.test_surface")),
            ("GUI mesh", ("-m", "unittest", "tools.oneflow_gui.test_mesh_validation")),
            ("GUI unrelated", ("-m", "unittest", "tools.oneflow_gui.test_parameter_model")),
            ("verification", ("-m", "unittest", "tools.verification.test_preflight")),
            ("verification discovery", ("-m", "unittest", "discover", "-s", "tools/verification", "-p", "test_*.py")),
            ("workflow discover", ("-m", "unittest", "discover", "-s", "tools", "-p", "test_workflow.py")),
            ("mesh catalog", ("tools/test_mesh_cases.py",)),
        )

    def test_stl_to_gmsh_selects_mesh_and_surface_dependents(self):
        selected = select_test_commands(("grid/stl_to_gmsh/converter.py",), self.registry)
        self.assertEqual(
            tuple(label for label, _ in selected),
            ("mesh geometry", "surface import", "GUI mesh"),
        )

    def test_verification_changes_select_verification_modules(self):
        selected = select_test_commands(("tools/verification/preflight.py",), self.registry)
        self.assertEqual(tuple(label for label, _ in selected), ("verification", "verification discovery"))

    def test_test_runner_registry_change_selects_registry_contract_tests(self):
        selected = select_test_commands(("tools/oneflow_gui/test_runner.py",), self.registry)
        self.assertEqual(
            tuple(label for label, _ in selected),
            ("verification",),
        )

    def test_gui_changes_select_gui_modules(self):
        selected = select_test_commands(("tools/oneflow_gui/home_window.py",), self.registry)
        self.assertEqual(tuple(label for label, _ in selected), ("GUI mesh", "GUI unrelated"))

    def test_unknown_code_area_falls_back_to_full_registry(self):
        selected = select_test_commands(("codes/solver/new_solver.cpp",), self.registry)
        self.assertEqual(tuple(label for label, _ in selected),
                         tuple(label for label, _ in self.registry))

    def test_docs_only_changes_select_no_python_tests(self):
        self.assertEqual(select_test_commands(("docs/README.md",), self.registry), ())

    def test_empty_change_set_selects_no_tests(self):
        self.assertEqual(select_test_commands((), self.registry), ())

    def test_mixed_known_and_unknown_changes_fall_back_to_full_suite(self):
        selected = select_test_commands(
            ("tools/verification/preflight.py", "new_component/runtime.cpp"),
            self.registry,
        )
        self.assertEqual(tuple(label for label, _ in selected),
                         tuple(label for label, _ in self.registry))


if __name__ == "__main__":
    unittest.main()
