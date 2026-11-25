# 90-degree circular pipe bend — geometry asset

This is the first curved-pipe geometry asset in the 3-D STL database. It is
a controlled analytic reference for the generic region-aware STL-to-Gmsh
volume-mesh route. Its generated MSH is registered as an existing-mesh input
for the Mesh-to-CFD Validation workflow.

## Canonical geometry

| Parameter | Default | Meaning |
|---|---:|---|
| Pipe radius (r) | 0.5 m | Circular cross-section radius |
| Bend centreline radius (R_c) | 1.0 m | Radius of the centreline arc |
| Bend angle (theta) | 90 deg | Direction change |
| Inlet straight length | 1.5 m | Straight upstream extension |
| Outlet straight length | 1.5 m | Straight downstream extension |
| Circumferential divisions | 48 | Points per cross-section ring |
| Bend divisions | 24 | Intervals along the circular arc |
| Inlet / outlet divisions | 12 / 12 | Axial intervals on straight sections |

The default dimensionless ratios are `Rc/D = 2` and
`Lin/D = Lout/D = 3`, where `D = 2r`.

The inlet centreline is parallel to +X. The centreline bends in the X-Y plane
and the outlet is tangent to the arc. The circular section is transported
with the centreline; it is not an extruded circle or a rotated straight STL.

For `0 <= phi <= theta`, the arc centreline is:

```text
C(phi) = (Rc sin(phi), Rc (1 - cos(phi)), 0)
T(phi) = (cos(phi), sin(phi), 0)
N(phi) = (-sin(phi), cos(phi), 0)
B    = (0, 0, 1)
```

The swept circular surface is:

```text
P(phi, alpha) = C(phi) + r cos(alpha) N(phi) + r sin(alpha) B
```

The inlet and outlet straight sections use the arc's endpoint tangent frames.
The theoretical enclosed fluid volume is:

```text
V = pi r^2 (Lin + Lout + Rc theta)
```

where `theta` is in radians. The triangulated STL volume is an approximation
and is checked against this analytic value.

## Boundary semantics

The shell has exactly three triangle patches:

- `wall`: swept lateral pipe wall;
- `inflow`: inlet disk that closes the fluid volume;
- `outflow`: outlet disk that closes the fluid volume.

The inlet and outlet disks are closure surfaces, not solid walls. The intended
Gmsh groups are dimension-2 `wall`, `inflow`, `outflow` and dimension-3
`fluid_domain`. The patch-to-triangle indices are stored in the Python
manifest. STL itself does not encode Gmsh Physical Groups. The current
mesh-project importer discovers the shared `region_manifest.json` beside the
STL, so keep one region-aware generated asset per output directory; the GUI's
default output directory already follows this rule.

## Generate outside the source tree

Python (binary STL plus checksum-bound JSON manifest with definition and generator hashes):

```powershell
python grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/python/generate_pipe_bend.py --output "$env:USERPROFILE/.oneflow/mesh-generation/generated_stl/pipe_bend_90_degree/pipe_bend_90_degree.stl"
```

C++17 on Windows (run from a Visual Studio Developer PowerShell):

```powershell
cl /std:c++17 /EHsc /O2 /Z7 /Fo"$env:USERPROFILE/.oneflow/mesh-generation/pipe_bend_generator.obj" /Fe"$env:USERPROFILE/.oneflow/mesh-generation/pipe_bend_generator.exe" grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/cpp/generate_pipe_bend.cpp
& "$env:USERPROFILE/.oneflow/mesh-generation/pipe_bend_generator.exe" --output "$env:USERPROFILE/.oneflow/mesh-generation/generated_stl/pipe_bend_90_degree/pipe_bend_90_degree_cpp.stl"
```

C++17 on Linux:

```bash
c++ -std=c++17 -O2 grid/stl_to_gmsh/3d/pipe_bend/90_degree/geometry/cpp/generate_pipe_bend.cpp -o "$HOME/.oneflow/mesh-generation/pipe_bend_generator"
"$HOME/.oneflow/mesh-generation/pipe_bend_generator" --output "$HOME/.oneflow/mesh-generation/generated_stl/pipe_bend_90_degree/pipe_bend_90_degree_cpp.stl"
```

Both generators accept the same parameter names: `--radius`,
`--bend-radius`, `--inlet-length`, `--outlet-length`,
`--angle-degrees`, `--circumferential-divisions`, `--bend-divisions`,
`--inlet-divisions`, `--outlet-divisions` and `--units` (`m`, `mm`, `cm`,
`in`). STL serialization differs
(binary Python, ASCII C++), so byte equality is not expected. Geometry,
triangle counts, boundary-patch counts and analytic dimensions are the
comparison contract.

## Mesh-to-CFD Validation registration

The generated-MSH profile is registered separately at
`../automatic_tetrahedra/`. It appears in the standalone 3-D Mesh-to-CFD
Validation tree after a project has produced its MSH. The profile discovers
completed projects by the `pipe_bend_90_degree` geometry identifier and
requires all four Physical Groups before conversion. This is intentionally an
existing-mesh consumer profile: geometry generation and volume meshing are
performed through the Geometry Library and STL-to-Gmsh workbench, while the
validation graph owns the MSH-to-CGNS-to-OneFLOW-to-CFD handoff.

## Verification

```powershell
python -m unittest grid.tests.python.stl_to_gmsh.test_pipe_bend_geometry
```

The tests check one connected watertight shell, manifold edge incidence,
consistent orientation, positive signed volume, analytic-volume agreement,
and complete non-overlapping patch assignment. The default Python mesh has
2,354 vertices and 4,704 triangles (2,304 wall, 48 inlet and 48 outlet
triangles).

## Scope boundary

The Python geometry generator writes three patch STL files and the existing
`oneflow.stl-surface-regions` `region_manifest.json` contract beside the closed
STL. `generate_volume_mesh.py` detects that manifest and dispatches to the
existing geometry-independent `_region_aware_closed_surface_volume_geo`
route. It merges the three patch files, coheres shared nodes, creates discrete
geometry, constructs one surface loop and one fluid volume, and assigns the
four required Physical Groups. Before meshing, the converter compares the
patch triangle geometry against the closed STL to verify an exact partition.
The adapter does not read the bend equations, angle or radius. Mesh-project
import copies the region manifest and all three patch files into the user's
project.

This is a generic STL-to-Gmsh route, not a pipe-bend-specific mesher. A real
Gmsh volume-mesh integration test is part of the existing Linux/Windows Gmsh
4.15.2 CI matrix. The generated-MSH profile is registered in Mesh-to-CFD
Validation and exercises CGNS and OneFLOW CFD smoke handoff after the user
creates the mesh in the GUI. Keep MSH 4.1 element/quality evidence with the
project; future work should promote the generation steps into a reproducible
mesh-producer Case only if the Case model can express the geometry parameters
and region-manifest handoff without duplicating the generic mesher. The C++17
geometry generator currently writes an ASCII STL and a parameter/patch-count
report; the region-manifest handoff is currently emitted by the Python
generator.
