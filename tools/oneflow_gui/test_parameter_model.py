"""Tests for the logical CFD/grid parameter model."""
from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.parameter_document import ParameterDocument
from tools.oneflow_gui.parameter_model import (
    ParameterApplicability,
    ParameterDomain,
    ParameterModel,
    classify_parameter,
)


class ParameterModelTest(unittest.TestCase):
    def _documents(self, root: Path, files: dict[str, str]):
        documents = []
        for name, content in files.items():
            path = root / name
            path.write_text(content, encoding="utf-8", newline="")
            documents.append(ParameterDocument.load(path))
        return documents

    def test_combines_shared_and_cfd_parameters_across_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "simu.txt": 'string simutask = "Solve";\n',
                "dimension.txt": "int dimension = 2;\n",
                "tolerence.txt": "Real ptTol = 1.0e-10;\n",
                "cfd.txt": (
                    "int maxSteps = 50;\n"
                    "Real global_dt = 0.001;\n"
                    "Real mach_ref = 0.729;\n"
                    "string visname = \"1eq-sa\";\n"
                ),
            })
            before = [document.original_bytes for document in documents]
            model = ParameterModel.from_documents(documents)

            self.assertEqual(model.task, "Solve")
            self.assertEqual(model.get("dimension").value, 2)
            self.assertEqual(model.get("ptTol").declared_type, "real")
            self.assertEqual(model.get("maxSteps").value, 50)
            self.assertEqual(model.get("global_dt").value, 0.001)
            self.assertEqual(model.get("visname").value, "1eq-sa")
            self.assertEqual(model.get("dimension").classification.domain, ParameterDomain.SHARED)
            self.assertEqual(model.get("dimension").classification.applicability, ParameterApplicability.BOTH)
            self.assertEqual(model.get("maxSteps").classification.domain, ParameterDomain.CFD)
            self.assertEqual(model.get("maxSteps").classification.applicability, ParameterApplicability.SOLVE)
            self.assertEqual(model.get("global_dt").classification.group, "Time stepping")
            self.assertEqual(
                {group.name for group in model.groups},
                {"General", "Run control", "Time stepping", "Flow and reference conditions", "Physical models"},
            )
            self.assertEqual([document.path.read_bytes() for document in documents], before)

    def test_grid_parameters_are_grouped_independently_of_file_names(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "simu.txt": 'string simutask = "Grid";\n',
                "dimension.txt": "int dimension = 2;\n",
                "grid.txt": (
                    "int gridObj = 1;\n"
                    'string sourceGridType = "plot3d";\n'
                    "Real gridScale = 1.0;\n"
                    'string sourceGridFileName = "grid/input.grd";\n'
                    "int parallelStrategy = 1;\n"
                ),
                "part.txt": (
                    "int npartproc = 4;\n"
                    "int partition_type = 1;\n"
                    'string ori_uns_file = "./grid/source.ofl";\n'
                ),
            })
            model = ParameterModel.from_documents(documents)

            self.assertEqual(model.task, "Grid")
            self.assertEqual(model.get("gridObj").classification.group, "Grid definition")
            self.assertEqual(model.get("sourceGridType").classification.group, "Grid formats and files")
            self.assertEqual(model.get("npartproc").classification.group, "Partitioning")
            self.assertEqual(model.get("dimension").classification.domain, ParameterDomain.SHARED)
            self.assertTrue(all(
                entry.classification.domain in {ParameterDomain.SHARED, ParameterDomain.GRID}
                for entry in model.entries
            ))


    def test_audited_cfd_case_settings_use_distinct_functional_groups(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "simu.txt": 'string simutask = "Solve";\\n',
                "cfd.txt": (
                    "int nEqu = 5;\\n"
                    "int nFieldSave = 50;\\n"
                    'string flowFile = "restart/flow.dat";\\n'
                    'string flowInputRoot = "default";\\n'
                    "int inflowType = 0;\\n"
                    "Real twall_dim = 300.0;\\n"
                    "Real gama_ref = 1.4;\\n"
                    "Real prl = 0.72;\\n"
                    "int lusgsFacePressureMode = 0;\\n"
                ),
            })
            model = ParameterModel.from_documents(documents)

            expected_groups = {
                "nEqu": "Solver system",
                "nFieldSave": "Output control",
                "flowFile": "Output files",
                "flowInputRoot": "Input and output files",
                "inflowType": "Boundary and inflow conditions",
                "twall_dim": "Boundary and inflow conditions",
                "gama_ref": "Thermodynamic reference state",
                "prl": "Transport properties",
                "lusgsFacePressureMode": "Numerical methods",
            }
            for key, expected_group in expected_groups.items():
                with self.subTest(parameter=key):
                    self.assertEqual(model.get(key).classification.group, expected_group)
                    self.assertEqual(
                        model.get(key).classification.applicability,
                        ParameterApplicability.SOLVE,
                    )


    def test_representative_case_manifests_have_no_unclassified_declarations(self):
        repository_root = Path(__file__).resolve().parents[2]
        case_paths = (
            Path("test/cfd/rae2822_roe_sa"),
            Path("test/grid/rae2822convert"),
        )
        for case_path in case_paths:
            with self.subTest(case=str(case_path)):
                script = repository_root / case_path / "script"
                control_text = (script / "control.txt").read_text(encoding="utf-8")
                documents = []
                for line in control_text.splitlines():
                    match = re.fullmatch(r'\s*"([^"]+)"\s*;\s*', line)
                    if match:
                        documents.append(ParameterDocument.load(script / match.group(1)))
                self.assertTrue(documents, f"No parameter documents found for {case_path}")
                model = ParameterModel.from_documents(documents)

                unclassified = [
                    entry.name for entry in model.entries
                    if entry.classification.domain == ParameterDomain.OTHER
                ]
                unclassified.extend(
                    declaration.name for declaration in model.unsupported_declarations
                    if classify_parameter(declaration.name).domain == ParameterDomain.OTHER
                )
                self.assertEqual(unclassified, [])

    def test_last_declaration_wins_and_shadowed_sources_remain_visible(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "first.txt": 'string simutask = "Solve";\nint maxSteps = 10;\n',
                "second.txt": "int maxSteps = 20;\n",
                "third.txt": "int MAXSTEPS = 30;\n",
            })
            model = ParameterModel.from_documents(documents)

            self.assertEqual(model.get("maxSteps").value, 30)
            self.assertEqual(model.get("maxSteps").source_file, "third.txt")
            self.assertEqual(len(model.shadowed_entries), 2)
            self.assertEqual(
                [entry.value for entry in model.shadowed_entries],
                [10, 20],
            )
            self.assertEqual(len([entry for entry in model.entries if entry.name.casefold() == "maxsteps"]), 1)

    def test_unknown_keys_are_retained_without_invented_classification(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "simu.txt": 'string simutask = "Solve";\n',
                "custom.txt": "Real vendorOption = 2.5;\n",
            })
            model = ParameterModel.from_documents(documents)
            entry = model.get("vendorOption")

            self.assertIsNotNone(entry)
            self.assertEqual(entry.value, 2.5)
            self.assertEqual(entry.declared_type, "real")
            self.assertEqual(entry.classification.domain, ParameterDomain.OTHER)
            self.assertEqual(entry.classification.group, "Unclassified")
            self.assertEqual(entry.classification.applicability, ParameterApplicability.UNKNOWN)

    def test_arrays_and_unrecognized_scalar_types_are_reported_as_unsupported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "grid.txt": "Real gridTrans[3] = {0.0, 0.0, 0.0};\ndouble vendorScale = 2.0;\n",
            })
            model = ParameterModel.from_documents(documents)

            self.assertEqual(
                [item.name for item in model.unsupported_declarations],
                ["gridTrans", "vendorScale"],
            )
            self.assertIn("Array declarations", model.unsupported_declarations[0].reason)
            self.assertIn("not supported", model.unsupported_declarations[1].reason)
            self.assertEqual(
                [item.line_number for item in model.unsupported_declarations],
                [1, 2],
            )

    def test_without_task_unknown_keys_stay_unclassified(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            documents = self._documents(root, {
                "custom.txt": "int vendorOption = 3;\n",
            })
            model = ParameterModel.from_documents(documents)
            entry = model.get("vendorOption")

            self.assertEqual(entry.classification.domain, ParameterDomain.OTHER)
            self.assertEqual(entry.classification.group, "Unclassified")
            self.assertEqual(entry.classification.applicability, ParameterApplicability.UNKNOWN)


if __name__ == "__main__":
    unittest.main()
