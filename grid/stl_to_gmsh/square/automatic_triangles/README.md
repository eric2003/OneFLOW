# Square with automatic TRI_3 refinement

## Purpose

This Case demonstrates the transition from a fixed two-triangle square to an
unstructured triangular mesh whose subdivision is selected by Gmsh. It is a
2D planar surface-meshing example. It does not hard-code a triangle count and
does not claim to support arbitrary curved or multiply connected STL surfaces.

## Input geometry

Both implementations use the same ASCII STL file, `input/model.stl`. It
contains two coplanar triangular facets forming the unit square. The diagonal
is part of the STL triangulation, not a requested CFD interior boundary.
The STL reader welds shared vertices and recovers the four-edge exterior
boundary.

The matching `script/layout.json` files define the surface Physical Group
`fluid_domain`, four boundary segments (`inflow`, `outflow`,
`symmetry`, and `farfield`), and `mesh_size = 0.25`. They deliberately
omit `transfinite_cells` and `recombine_surface`.

## Construction and implementation

### Native C++

The Case is launched through the OneFLOW Grid task using
`script/dimension.txt`, `script/simu.txt`, and `script/grid.txt`.
The grid parameters identify the STL input, the JSON layout, the Gmsh
executable, and `results/square_automatic_triangles.msh` as the artifact.
The native conversion path recovers the planar boundary, emits a Gmsh
geometry script, sets the target characteristic length, invokes Gmsh, and
reads the MSH 4.1 result for topology and boundary validation.

### Python reference

The Python Case uses the independent
`grid/stl_to_gmsh/converter.py` implementation through
`grid/stl_to_gmsh/run_python_case.py`. Its `expected.json` describes
stable mesh invariants and the shared artifact name. It does not call the
native C++ producer.

## How refinement is controlled

Gmsh receives `Mesh.MeshSizeMin` and `Mesh.MeshSizeMax` equal to the
configured `mesh_size`, with point and boundary size propagation enabled.
Reducing `mesh_size` requests a finer mesh. The exact node and triangle
counts remain mesher outputs and may vary with geometry, algorithm, and Gmsh
version. The progressive-refinement integration test compares several
mesh sizes and checks the expected trend without prescribing exact counts.

This is intentionally different from `two_triangles`, whose transfinite
configuration fixes a topology and an exact two-cell result.

## Artifact and acceptance contract

The canonical artifact is
`results/square_automatic_triangles.msh` for both implementations.
Validation checks:

- more than four nodes and at least four exterior boundary LINE_2 elements;
- one or more TRI_3 cells and no QUAD_4 cells;
- only the expected line and triangle element families;
- boundary LINE_2 edges agree with the boundary edges of the triangle mesh;
- valid edge incidence and the Euler characteristic of a simply connected
  planar square;
- all five Physical Groups are retained.

The Euler check is specific to this simply connected square fixture; it is
not a general proof for domains with holes or multiple connected components.

## Workflow and evidence

The Case is registered in `grid/case_catalog.json` and is selectable under
2D → STL to Gmsh in Run Development Cases, Python reference cases, and
Mesh-to-CFD Validation. The latter follows the shared artifact-bound
Gmsh → CGNS → OneFLOW handoff. For this simple square, the downstream CFD
consumer is a short runtime smoke test only; it is not a numerical accuracy
benchmark.

Run the isolated producer Cases from the repository root:

```powershell
python grid/stl_to_gmsh/develop_test/run_cases.py --oneflow-executable build/bin/Debug/OneFLOW.exe --gmsh gmsh
```

The development runner stages inputs under `%USERPROFILE%/.oneflow` on
Windows or `$HOME/.oneflow` on Linux. Never run the checked-in Case directory
in-place. The automated tests live under `tests/python/grid` and
`tests/geometry`; they do not create generated files in the source Case.

## Scope limits

This Case currently targets one connected, axis-aligned rectangular XY
boundary and a simply connected planar surface. It does not establish support
for arbitrary non-planar STL, multiple holes, self-intersections, curved CAD
entities, or general 3D volume meshing. The 3D cube prototype is a separate
Case family member with a separate volume topology contract.
