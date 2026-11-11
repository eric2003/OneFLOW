# Square with Two TRI_3 Cells

This pair of Cases uses the same unit-square STL as the single-quad square
Case. The STL boundary is unchanged; only the surface recombination policy
changes.

- Geometry: unit square in the XY plane, x and y from 0 to 1.
- Transfinite boundary: one segment on each of the four edges.
- Surface recombination: disabled.
- Surface Physical Group: fluid_domain.
- Boundary Physical Groups: inflow, outflow, symmetry and farfield.

Expected MSH 4.1 output:

| Quantity | Expected |
|---|---:|
| Nodes | 4 |
| Boundary LINE_2 elements | 4 |
| TRI_3 cells | 2 |
| QUAD_4 cells | 0 |

The diagonal between the two generated triangles is an internal cell interface,
not a boundary Physical Group. The four exterior Physical Groups remain
identical to the single-quad Case.

The native C++ Case uses the standard OneFLOW Grid task. The Python reference
uses the independent STL-to-Gmsh implementation. Both are staged outside the
repository before execution; generated outputs do not belong in this directory.
