# STL to Gmsh: Square 2 x 2

This native C++ OneFLOW Grid Case uses the exact same STL geometry as
../square/cpp/input/model.stl. The two STL files are byte-identical. Only the
transfinite mesh layout changes.

## Geometry and mesh

- Geometry: unit square in the XY plane, x and y from 0 to 1.
- Surface Physical Group: fluid_domain.
- Left boundary: inflow.
- Right boundary: outflow.
- Bottom boundary: symmetry.
- Top boundary: farfield.
- Transfinite cells: 2 in x and 2 in y.

Expected MSH 4.1 output:

| Quantity | Expected |
|---|---:|
| Nodes | 9 |
| Boundary LINE_2 elements | 8 |
| QUAD_4 cells | 4 |
| Triangle elements | 0 |

The four boundary Physical Groups each contain two line elements. The STL
triangulation diagonal is not a CFD cell interface.

## Parameters

- script/simu.txt selects the Grid task.
- script/dimension.txt selects two dimensions.
- script/grid.txt selects STL -> Gmsh and the output path.
- script/layout.json defines the transfinite cell counts and physical edges.

The Case is staged outside the repository by the development runner. Generated
MSH files and reports do not belong in this source directory.
