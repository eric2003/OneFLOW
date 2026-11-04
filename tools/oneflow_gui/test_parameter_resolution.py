"""Tests for file-independent Case parameter override resolution."""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.parameter_document import ParameterDocument
from tools.oneflow_gui.parameter_resolution import (
    ParameterOverride,
    resolve_parameter_overrides,
)


class ParameterResolutionTest(unittest.TestCase):
    def _documents(self, root: Path, files: dict[str, str]):
        result = []
        for name, content in files.items():
            path = root / name
            path.write_text(content, encoding="utf-8", newline="")
            result.append(ParameterDocument.load(path))
        return result

    def test_max_steps_override_is_independent_of_legacy_filename(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "control.txt": 'string simutask = "Solve";\n',
                "run-control.txt": (
                    "int maxSteps = 50; // Total simulation steps\n"
                    "int iexitflag = 0;\n"
                ),
            })
            original = [item.path.read_bytes() for item in documents]

            plan = resolve_parameter_overrides(
                documents,
                [ParameterOverride("maxSteps", 35000, "production-profile")],
            )

            self.assertEqual(len(plan.changes), 1)
            change = plan.changes[0]
            self.assertEqual(change.previous_value, 50)
            self.assertEqual(change.requested_value, 35000)
            self.assertEqual(Path(change.source_path).name, "run-control.txt")
            self.assertIn("maxSteps: 50 -> 35000", plan.preview()[0])
            self.assertEqual([item.path.read_bytes() for item in documents], original)

            plan.apply()
            self.assertIn(
                "int maxSteps = 35000; // Total simulation steps",
                documents[1].path.read_text(encoding="utf-8"),
            )
            self.assertEqual(
                documents[0].path.read_text(encoding="utf-8"),
                'string simutask = "Solve";\n',
            )

    def test_effective_duplicate_declaration_is_updated_in_owning_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "first.txt": "int maxSteps = 10;\n",
                "later.txt": "int maxSteps = 50;\n",
            })
            plan = resolve_parameter_overrides(documents, {"MAXSTEPS": 35000})
            self.assertEqual(Path(plan.changes[0].source_path).name, "later.txt")
            plan.apply()
            self.assertEqual(documents[0].get("maxSteps").value, 10)
            self.assertEqual(documents[1].get("maxSteps").value, 35000)

    def test_missing_key_is_rejected_without_modifying_any_document(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {"run.txt": "int maxSteps = 50;\n"})
            before = [item.path.read_bytes() for item in documents]
            with self.assertRaises(KeyError):
                resolve_parameter_overrides(documents, {"maxTime": 0.2})
            self.assertEqual([item.path.read_bytes() for item in documents], before)

    def test_duplicate_override_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {"run.txt": "int maxSteps = 50;\n"})
            with self.assertRaises(ValueError):
                resolve_parameter_overrides(documents, [
                    ParameterOverride("maxSteps", 100),
                    ParameterOverride("MAXSTEPS", 200),
                ])

    def test_multi_file_apply_updates_each_owning_document(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "run.txt": "int maxSteps = 50;\\n",
                "time.txt": "Real global_dt = 0.001;\\n",
            })
            plan = resolve_parameter_overrides(
                documents, {"maxSteps": 35000, "global_dt": 0.002}
            )

            plan.apply()

            self.assertEqual(documents[0].get("maxSteps").value, 35000)
            self.assertEqual(documents[1].get("global_dt").value, 0.002)
            self.assertEqual(
                documents[0].path.read_text(encoding="utf-8"),
                "int maxSteps = 35000;\\n",
            )
            self.assertEqual(
                documents[1].path.read_text(encoding="utf-8"),
                "Real global_dt = 0.002;\\n",
            )

    def test_multi_file_write_failure_rolls_back_prior_replacements(self):
        from unittest.mock import patch

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "run.txt": "int maxSteps = 50;\\n",
                "time.txt": "Real global_dt = 0.001;\\n",
            })
            plan = resolve_parameter_overrides(
                documents, {"maxSteps": 35000, "global_dt": 0.002}
            )
            original = [document.path.read_bytes() for document in documents]

            with patch(
                "tools.oneflow_gui.parameter_resolution.os.replace",
                side_effect=[None, OSError("simulated second-file replacement failure"), None],
            ):
                with self.assertRaisesRegex(RuntimeError, "rolled back"):
                    plan.apply()

            self.assertEqual(
                [document.path.read_bytes() for document in documents], original
            )
            self.assertEqual(
                [document.original_bytes for document in documents], original
            )
            self.assertEqual(list(root.glob(".*.oneflow-tmp")), [])

    def test_external_change_is_detected_before_any_file_is_written(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "a.txt": "int maxSteps = 50;\n",
                "b.txt": "Real global_dt = 0.001;\n",
            })
            plan = resolve_parameter_overrides(
                documents, {"maxSteps": 35000, "global_dt": 0.002}
            )
            documents[1].path.write_text(
                "Real global_dt = 0.003;\n", encoding="utf-8"
            )
            before_first = documents[0].path.read_bytes()
            with self.assertRaises(RuntimeError):
                plan.apply()
            self.assertEqual(documents[0].path.read_bytes(), before_first)


if __name__ == "__main__":
    unittest.main()
