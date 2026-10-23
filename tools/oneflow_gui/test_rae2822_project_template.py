"""Tests for the RAE2822 multi-case project template."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.oneflow_gui.project_templates import create_project


class Rae2822ProjectTemplateTest(unittest.TestCase):
    def test_copies_cases_and_dependency_workflow_to_workspace(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            repository = root / "repository"
            grid_case = repository / "grid" / "mesh_workflow"
            flow_case = repository / "test" / "cfd" / "rae2822_roe_sa"
            example = repository / "tools" / "oneflow_gui" / "examples"
            for directory in (grid_case / "script", flow_case / "script", example):
                directory.mkdir(parents=True)
            (grid_case / "script" / "grid.txt").write_text(
                'string targetGridFileName = "results/rae2822_2d.ofl";\n',
                encoding="utf-8",
            )
            (flow_case / "script" / "cfd.txt").write_text(
                'string gridFileName = "grid/rae2822_2d.ofl";\n',
                encoding="utf-8",
            )
            workflow = {
                "schema_version": 1,
                "workflow_id": "rae2822-grid-to-solver",
                "name": "RAE2822 grid conversion and CFD",
                "steps": [{
                    "type": "case_graph",
                    "cases": [
                        {
                            "id": "rae2822grid",
                            "directory": "rae2822grid",
                            "outputs": [{
                                "id": "rae2822_grid",
                                "path": "results/rae2822_2d.ofl",
                                "parameter_file": "script/grid.txt",
                                "setting": "targetGridFileName",
                            }],
                        },
                        {
                            "id": "rae2822_roe_sa",
                            "directory": "rae2822_roe_sa",
                            "depends_on": ["rae2822grid"],
                            "inputs": [{
                                "artifact": "rae2822_grid",
                                "setting": "gridFileName",
                                "parameter_file": "script/cfd.txt",
                            }],
                        },
                    ],
                }],
            }
            (example / "rae2822_multicase.json").write_text(
                json.dumps(workflow), encoding="utf-8"
            )

            project = create_project(
                repository, "rae2822_multicase", "RAE2822Workflow", root / "workspace"
            )

            self.assertTrue((project.root / "rae2822grid" / "script" / "grid.txt").is_file())
            self.assertTrue((project.root / "rae2822_roe_sa" / "script" / "cfd.txt").is_file())
            self.assertTrue((project.root / ".oneflow-project.json").is_file())
            self.assertTrue((project.root / "README.md").is_file())
            copied = json.loads((project.root / "workflow.json").read_text(encoding="utf-8"))
            cases = copied["steps"][0]["cases"]
            self.assertEqual(cases[1]["depends_on"], ["rae2822grid"])
            self.assertEqual(cases[0]["outputs"][0]["path"], "results/rae2822_2d.ofl")
            self.assertEqual(cases[1]["inputs"][0]["artifact"], "rae2822_grid")
            self.assertTrue((grid_case / "script" / "grid.txt").is_file())


if __name__ == "__main__":
    unittest.main()
