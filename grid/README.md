# Grid Cases

Each directory under `grid/` is a concrete mesh project or case. Case files
describe inputs, parameters, expected results, and how to run the workflow.
Generated grids, reports, and visualization artifacts belong in an external
workspace, not in the checked-in case definition.

## Case catalog

`case_catalog.json` is the machine-readable discovery index for tools and a
future GUI. It gives each case a stable ID and locates its runtime entry point,
layout, and optional richer profile. A GUI should discover cases from this
catalog instead of hard-coding case names or directory paths.

The `minimal_rectangle`, `minimal_cube`, and `two_connected_cubes` cases now
have versioned `case.json` profiles. Each profile describes the generation
dimension, layout, runtime parameter-file list, inspection outputs, and
expected mesh counts and bounds. The catalog remains an index; profile
validation checks that its layout agrees with the catalog and that its listed
parameter files exist.

Use the standard-library command-line tool to discover and validate cases:

```sh
python tools/mesh_cases.py list
python tools/mesh_cases.py show minimal_cube
python tools/mesh_cases.py validate
python tools/mesh_cases.py validate two_connected_cubes
python tools/test_mesh_cases.py
```

The catalog and profiles are not yet the authoritative store for every runtime
parameter. Existing `script/*.txt` files remain authoritative for OneFLOW, so
standalone execution and older workflows are unchanged. The next migration
stage is to generate selected legacy parameter files from explicitly modeled
profile settings and compare generated values with the existing inputs. Values
not yet modeled must be preserved rather than silently discarded. Do not
duplicate parameter values in the catalog.

## Cases

| Case | Purpose | Configuration and execution |
|---|---|---|
| `minimal_rectangle` | Generate a small structured grid, inspect geometry and topology, and export SVG/VTU views | `case.json`, `run.ps1` |
| `minimal_cube` | Generate one 3D block from a readable point-line-face-block layout | `case.json`, `script/control.txt` |
| `two_connected_cubes` | Generate two adjacent 3D blocks with one shared interface face | `case.json`, `script/control.txt` |
| `mesh_workflow` | Convert the RAE 2822 Plot3D grid to OneFLOW format and write a metadata summary | `script/control.txt`, parameter files, `tools/validate_case.py` |

## Target integration model

- Keep the case profile independent of GUI widgets and OneFLOW runtime globals.
- Let GUI, CLI, and tests use the same case-profile reader and workflow service.
- Generate legacy parameter files only at the case/workspace boundary; do not
  make OneFLOW depend on a GUI or the case catalog.
- Launch the normal OneFLOW executable for actual generation and conversion.
- Validate produced geometry and topology against case expectations, not just
  process exit codes or file existence.
- Keep generated grids, reports, and visualization artifacts in an external
  workspace.

## Organization rules

- Keep case-specific parameters and expected results with the case.
- Share code only after multiple concrete cases demonstrate the same need.
- Preserve existing parameter-file compatibility during migration.

## Desktop GUI

The optional PySide6 desktop runner discovers cases from the same catalog,
runs selected cases sequentially, streams process logs, and stores each run in
an external workspace so source case files are not overwritten.

From the repository root:

```powershell
python -m pip install -r tools/requirements-gui.txt
python -m tools.oneflow_gui
python -m unittest tools.oneflow_gui.test_services
```

See `tools/oneflow_gui/README.md` for module boundaries and run behavior.
The GUI currently reports process status; geometry and topology validation
remain a separate workflow step and are not inferred from exit code alone.
