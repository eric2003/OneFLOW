# STL to Gmsh - 2-D structured mesh generation

This is the first runnable STL-to-Gmsh mesh-generation workflow. It reads a
planar square STL boundary, recovers its four outer edges, creates a Gmsh
geometry, generates an ASCII MSH 4.1 mesh, and validates the mesh topology and
Physical Groups.

The single-cell STL contains two triangles that describe the unit-square
boundary. The two-cell STL describes a 2 x 1 rectangle with two unit-square
cells. STL diagonals are only triangulation edges; they are not interpreted
as CFD boundaries or as requested cell interfaces. The requested cell layout
is explicit in the sidecar files under `script/`.

## Examples

| Example | Structured layout | Nodes | QUAD_4 cells | LINE_2 boundary elements |
|---|---:|---:|---:|---:|
| `single_quad` | 1 x 1 | 4 | 1 | 4 |
| `two_quads` | 2 adjacent 1 x 1 squares | 6 | 2 | 6 |

The two-cell example uses a 2 x 1 rectangle split at x = 1, so each cell is
itself a 1 x 1 square. It has one internal interface shared by the two
quadrilaterals; that interface is not an external boundary condition.

## Boundary contract

The four outer edges are explicitly assigned in the JSON sidecar:

| Edge | Gmsh Physical Curve | Intended boundary |
|---|---|---|
| Left | `inflow` | Inflow |
| Right | `outflow` | Outflow |
| Bottom | `symmetry` | Symmetry |
| Top | `farfield` | Farfield |

The 2-D region is assigned to the `fluid_domain` Physical Surface. These
names are intentionally lowercase because the existing Gmsh-to-CGNS adapter
maps `inflow`, `outflow`, `symmetry`, and `farfield` to the corresponding
CGNS BC types. The names identify boundary types; they do not define
solver-specific boundary state values.

## Run

Gmsh must be installed and its executable available on PATH. The workflow
writes only to the staged case directory, never to the repository source tree.

PowerShell:

```powershell
$case = Join-Path $env:USERPROFILE ".oneflow/cases/stl_to_gmsh"
New-Item -ItemType Directory -Force (Split-Path $case) | Out-Null
if (Test-Path $case) { Remove-Item -Recurse -Force $case }
Copy-Item -Recurse -Force ".\grid\stl\stl_to_gmsh" $case
python ".\grid\stl\stl_to_gmsh\run_case.py" --case-directory $case
```

Linux:

```sh
case_dir="$HOME/.oneflow/cases/stl_to_gmsh"
mkdir -p "$(dirname "$case_dir")"
rm -rf "$case_dir"
cp -R grid/stl/stl_to_gmsh "$case_dir"
python grid/stl/stl_to_gmsh/run_case.py --case-directory "$case_dir"
```

To generate only the two-cell example, add `--example two_quads`. To use a
non-default Gmsh installation, pass `--gmsh /path/to/gmsh` (or the Windows
executable path).

## Outputs

The staged case receives:

- `results/single_quad.geo` and `results/single_quad.msh`
- `results/two_quads.geo` and `results/two_quads.msh`
- One JSON report per example, including node/element counts and Physical
  Group names.

The generated MSH files are Gmsh 4.1 ASCII, matching the current OneFLOW
Gmsh-to-CGNS reader. They can be used as input to the existing
`gmsh -> cgns` conversion route. The next workflow step is to hand these
generated MSH files to the standard OneFLOW CGNS/mesh conversion Case.

## Current scope

This first generator supports one connected, closed, axis-aligned planar
quadrilateral STL boundary and structured quadrilateral layouts. The
`transfinite_cells.x` and `transfinite_cells.y` values define the number of
cells in each direction. Arbitrary polygonal boundaries, holes, multiple
regions, curved boundaries, and 3-D surface/volume meshing are separate
capabilities and are not inferred from STL triangle edges.
