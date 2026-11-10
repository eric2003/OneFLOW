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

The foundational 2D and 3D layout cases have versioned `case.json` profiles. Each profile describes the generation dimension, layout, runtime parameter-file list, inspection outputs, and expected mesh counts and bounds. The flat-plate cases also provide profiles
with mesh expectations; `turbulent_flat_plate` follows the laminar reference
layout but uses 81 points in the wall-normal direction. The catalog remains an
index; profile validation checks that its layout agrees with the catalog and
that its listed parameter files exist.

Use the standard-library command-line tool to discover and validate cases:

```sh
python tools/mesh_cases.py list
python tools/mesh_cases.py show minimal_cube
python tools/mesh_cases.py validate
python tools/mesh_cases.py validate two_blocks_shared_face_corners_only
python tools/mesh_cases.py validate two_blocks_shared_edge_endpoints_only
python tools/mesh_cases.py validate two_blocks_shared_edge_one_interior_node
python tools/mesh_cases.py validate turbulent_flat_plate
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
| `minimal_rectangle` | Generate a small structured grid, inspect geometry and topology, and export SVG/VTU views | `grid/mesh-layout/2d/minimal_rectangle` |
| `minimal_rectangle_vertical_midpoints` | Rectangle layout with one midpoint on each vertical side | `grid/mesh-layout/2d/minimal_rectangle_vertical_midpoints` |
| `minimal_cube` | Generate one 3D block from a readable point-line-face-block layout | `grid/mesh-layout/3d/minimal_cube` |
| `two_blocks_shared_face_corners_only` | Two 3D blocks sharing a face represented by its four corner nodes | `grid/mesh-layout/3d/two_blocks_shared_face_corners_only` |
| `two_blocks_shared_edge_endpoints_only` | Minimum 2D two-block case: one QUAD per block and two shared-edge endpoints | `grid/mesh-layout/2d/two_blocks_shared_edge_endpoints_only` |
| `two_blocks_shared_edge_one_interior_node` | Two 2D blocks whose shared edge contains one interior node | `grid/mesh-layout/2d/two_blocks_shared_edge_one_interior_node` |
| `mesh_workflow` | Convert the RAE 2822 Plot3D grid to OneFLOW format and write a metadata summary | `script/control.txt`, parameter files, `tools/validate_case.py` |
| `gmsh_to_cgns` | Manual GUI example for Gmsh MSH 4.1 ASCII to CGNS conversion | `script/control.txt`, `grid/square.msh`, `README.md` |
| `gmsh_to_oneflow` | Manual GUI example for Gmsh MSH 4.1 ASCII through CGNS to a native OneFLOW grid | `script/control.txt`, `grid/single_quad.msh`, `README.md` |
| `two_quads` | Self-contained two-cell Gmsh conversion example with an explicitly ungrouped internal interface | `grid/conversion/cases/two_quads` |
| `stl_to_gmsh_square`, `stl_to_gmsh_two_square_cells` | Native C++ OneFLOW STL-to-Gmsh Cases with independent Python references | `grid/stl_to_gmsh/<case>/cpp`, `grid/stl_to_gmsh/<case>/python` |
| `laminar_flat_plate` | Reference two-block flat-plate layout template with 85 wall-normal points | `grid/mesh-layout/2d/laminar_flat_plate` |
| `laminar_flat_plate_10_edges` | Two-block laminar flat-plate layout with ten geometric line segments and 41 points across the first 0.02 near-wall region | `grid/mesh-layout/2d/laminar_flat_plate_10_edges` |
| `turbulent_flat_plate` | Flat-plate layout variant with 81 wall-normal points and updated expected topology counts | `grid/mesh-layout/2d/turbulent_flat_plate` |

The native `stl_to_gmsh` Cases run through the OneFLOW Grid task and use the configured Gmsh CLI to generate MSH 4.1. Their outputs can feed the existing Gmsh-to-CGNS/OneFLOW conversion paths. Each problem keeps its Native C++ and Python Reference implementations side by side under `grid/stl_to_gmsh/<case>/cpp` and `grid/stl_to_gmsh/<case>/python`.

## Mesh layout templates

Foundational layout cases are grouped by dimension:

- `2d/`: rectangle baselines, the two-block shared-edge progression, and flat-plate layouts.
- `3d/`: the minimal cube and the two-block shared-face baseline.

The two-block 2D examples are deliberately separate. The endpoint-only case
has two nodes on the shared edge and one QUAD per block. The next case adds
one interior node to that edge. Each case tests only its stated topology
condition; more complex geometry and boundary-condition behavior belong in
separate examples.

Each template remains a self-contained case with its own profile where
applicable, runtime parameter files, layout source, and editable specification.
Generated artifacts belong in an external workspace.

## Conversion Case organization

Conversion examples are grouped under grid/conversion/cases/<case-id>/.
The directory name describes the mesh/topology, not a source/target format
pair. Each runnable Case owns its source mesh and runtime configuration so it
can be copied and executed independently. The central
grid/interoperability/conversion_matrix.json records implemented routes
and the concrete Cases that cover them. Existing route-named Cases remain
registered for backward compatibility; new examples use topology-oriented
names.

The first example, two_quads, contains two adjacent QUAD_4 cells, six
nodes, six exterior LINE_2 elements, and one internal interface line. Its
four exterior Physical Groups demonstrate that a single named group can span
multiple geometric entities. See grid/conversion/README.md and the
Grid conversion Cases documentation for the current collection and scope.

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

## Runtime input and output roots

Runtime path roots are optional settings in the existing `script/*.txt`
parameter files. Existing cases that omit them continue to use the active case
directory.

The general settings are `projectRoot`, `inputRoot`, and `outputRoot`.
Mesh and flow data can then override them independently with
`meshInputRoot`, `flowInputRoot`, `meshOutputRoot`, and `flowOutputRoot`.
A missing domain-specific setting inherits its corresponding general root.
An empty setting or `default` selects that inherited value; `case` selects
the active case; `project` selects `projectRoot`; absolute paths are used
directly; other relative IO roots are resolved from the active case. A relative
`projectRoot` is resolved from the process working directory.

For example, a conversion case can read an initial grid from a shared mesh
directory and write the converted grid to its own output directory, while a
CFD case can read that converted grid and write flow results elsewhere:

```text
projectRoot = ../OneFLOWProject
meshInputRoot = ../shared-mesh-case
meshOutputRoot = ../mesh-output-case
flowInputRoot = ../mesh-output-case
flowOutputRoot = ../flow-output-case
```

The root settings do not replace existing file-name settings such as
`sourceGridFileName` or `targetGridFileName`. A root is a base directory, not
a replacement for the `grid/` or `results/` subdirectory already present in a
file-name value. Relative file names are appended to the applicable root,
while absolute file paths remain absolute. Keep existing file-name values
unchanged for compatibility; when moving a case to a different artifact root,
check that the configured relative file names still match the producer output.
Each run writes an `io_paths.txt` diagnostic beneath each configured output
root (only once when both roots resolve to the same path). When both output
roots are the case root, the legacy-friendly location `results/io_paths.txt` is
used. The report records the effective roots, so they can be inspected without
inferring them from parameters.

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
