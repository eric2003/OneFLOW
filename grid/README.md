# Grid Cases

Each directory under `grid/` is a concrete mesh project or case. Case files
describe inputs, parameters, expected results, and how to run the workflow.
Generated grids, reports, and visualization artifacts belong in an external
workspace, not in the checked-in case definition.

## Cases

| Case | Purpose | Configuration and execution |
|---|---|---|
| `minimal_rectangle` | Generate a small structured grid, inspect geometry and topology, and export SVG/VTU views | `case.json`, `run.ps1` |
| `mesh_workflow` | Convert the RAE 2822 Plot3D grid to OneFLOW format and write a metadata summary | `script/control.txt`, parameter files, `tools/validate_case.py` |

## Organization rules

- Keep case-specific parameters and expected results with the case.
- Run calculations from a separate workspace and direct all generated output there.
- Share code only after multiple concrete cases demonstrate the same need.
- Keep the validation tied to the case's real output, not only to process exit
  codes or the existence of files.
