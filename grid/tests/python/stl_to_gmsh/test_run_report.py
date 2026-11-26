from __future__ import annotations
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock
from grid.stl_to_gmsh import generate_volume_mesh
from grid.stl_to_gmsh.run_report import build_run_report, format_run_summary, stage_result, write_run_report

class RunReportTest(unittest.TestCase):
    def test_report_rejects_unknown_stage_status(self):
        with self.assertRaisesRegex(ValueError, "Unsupported run-stage status"):
            stage_result("MAYBE")

    def test_summary_distinguishes_not_run_and_not_checked(self):
        report = build_run_report("STL -> Gmsh", {
            "geometry": stage_result("PASS"),
            "gmsh": stage_result("FAIL", "Gmsh exited with code 1"),
            "mesh_structure": stage_result("NOT_RUN"),
            "mesh_quality": stage_result("NOT_CHECKED", "No acceptance check is implemented"),
        }, 2)
        summary = format_run_summary(report)
        self.assertIn("STL -> Gmsh: FAILED", summary)
        self.assertIn("[PASS] geometry", summary)
        self.assertIn("[FAIL] gmsh", summary)
        self.assertIn("[NOT RUN] mesh_structure", summary)
        self.assertIn("[NOT CHECKED] mesh_quality", summary)

    def test_summary_keeps_file_inventory_in_json_not_console_output(self):
        report = build_run_report(
            "Mesh-to-CFD Validation",
            {"mesh_generation": stage_result("PASS")},
            0,
            artifacts={
                "run_directory": "run",
                "files": [{
                    "path": "stl_to_gmsh/results/mesh.msh",
                    "sha256": "a" * 64,
                    "size_bytes": 1024,
                }],
            },
        )
        summary = format_run_summary(report)
        self.assertIn("run_directory: run", summary)
        self.assertNotIn("  files:", summary)
        self.assertNotIn("sha256", summary)

    def test_success_report_cannot_contain_a_failed_stage(self):
        with self.assertRaisesRegex(ValueError, "cannot contain a failed stage"):
            build_run_report("STL -> Gmsh", {"gmsh": stage_result("FAIL")}, 0)

    def test_report_is_written_as_valid_json(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "run.report.json"
            report = build_run_report("STL -> Gmsh", {"input": stage_result("PASS")}, 0)
            write_run_report(path, report)
            self.assertEqual(json.loads(path.read_text(encoding="utf-8")), report)
            self.assertEqual(list(Path(temporary).glob("*.tmp")), [])

class GenerateVolumeMeshReportTest(unittest.TestCase):
    def _mesh_summary(self):
        return {"format":"4.1","nodes":10,"elements":20,"triangle_elements":8,
                "tetrahedron_elements":12,
                "physical_groups":[{"dimension":2,"name":"wall","element_count":8},{"dimension":3,"name":"fluid_domain","element_count":12}],
                "triangle_quality":{"minimum_quality":0.8}}

    def test_success_prints_short_summary_and_writes_stage_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); stl=root/"input.stl"; stl.write_text("test",encoding="utf-8")
            output=root/"results"/"mesh.msh"
            with (mock.patch.object(generate_volume_mesh,"convert",return_value={"triangles":4}),
                  mock.patch.object(generate_volume_mesh.subprocess,"run"),
                  mock.patch.object(generate_volume_mesh,"inspect_msh4",return_value=self._mesh_summary()),
                  mock.patch("builtins.print") as printed):
                report=generate_volume_mesh.generate(stl,output,"gmsh",0.1,0.2)
            saved=json.loads(output.with_suffix(".report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["status"],"success")
            self.assertEqual(saved["exit_code"],0)
            self.assertTrue(saved["validated"])
            for stage in ("geometry","gmsh","mesh_structure","physical_groups"):
                self.assertEqual(saved["stages"][stage]["status"],"PASS")
            self.assertEqual(saved["stages"]["mesh_quality"]["status"],"NOT_CHECKED")
            self.assertEqual(saved["stages"]["cfd_validation"]["status"],"NOT_RUN")
            rendered="\n".join(str(c.args[0]) for c in printed.call_args_list)
            self.assertIn("STL -> Gmsh: SUCCESS",rendered)
            self.assertIn("[NOT CHECKED] mesh_quality",rendered)
            self.assertNotIn('"schema"',rendered)

    def test_missing_gmsh_writes_preflight_failure_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); stl=root/"input.stl"; stl.write_text("test",encoding="utf-8")
            output=root/"results"/"mesh.msh"
            with mock.patch("builtins.print") as printed:
                with self.assertRaisesRegex(ValueError, "Gmsh was not found"):
                    generate_volume_mesh.generate(stl,output,None,0.1,0.2)
            saved=json.loads(output.with_suffix(".report.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["status"],"failed")
            self.assertEqual(saved["stages"]["input"]["status"],"FAIL")
            self.assertEqual(saved["stages"]["geometry"]["status"],"NOT_RUN")
            self.assertIn("Gmsh was not found",saved["error"]["message"])
            self.assertIn("STL -> Gmsh: FAILED",printed.call_args.args[0])

    def test_gmsh_failure_writes_failure_report_and_marks_later_stages_not_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); stl=root/"input.stl"; stl.write_text("test",encoding="utf-8")
            output=root/"results"/"mesh.msh"
            failure=subprocess.CalledProcessError(1,["gmsh"],stderr="Error: meshing failed\n")
            with (mock.patch.object(generate_volume_mesh,"convert",return_value={"triangles":4}),
                  mock.patch.object(generate_volume_mesh.subprocess,"run",side_effect=failure),
                  mock.patch("builtins.print") as printed):
                with self.assertRaises(subprocess.CalledProcessError):
                    generate_volume_mesh.generate(stl,output,"gmsh",0.1,0.2)
            saved=json.loads(output.with_suffix(".report.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["status"],"failed")
            self.assertEqual(saved["exit_code"],2)
            self.assertFalse(saved["validated"])
            self.assertEqual(saved["stages"]["geometry"]["status"],"PASS")
            self.assertEqual(saved["stages"]["gmsh"]["status"],"FAIL")
            self.assertEqual(saved["stages"]["mesh_structure"]["status"],"NOT_RUN")
            self.assertIn("meshing failed",saved["error"]["message"])
            rendered="\n".join(str(c.args[0]) for c in printed.call_args_list)
            self.assertIn("STL -> Gmsh: FAILED",rendered)
            self.assertIn("[FAIL] gmsh",rendered)

if __name__ == "__main__":
    unittest.main()
