# Square STL mesh variants

The square geometry is shared by the single-quad and two-triangle Cases.
The STL is the same two-facet representation in both Cases; its diagonal is
removed when extracting the outer boundary and is not itself a CFD cell
interface.

| Variant | Layout | Surface recombination | Output |
|---|---|---|---|
| square/single_quad/cpp and square/single_quad/python | 1 x 1 | enabled | 1 QUAD_4 |
| two_triangles/cpp and two_triangles/python | 1 x 1 | disabled | 2 TRI_3 |

Both variants have four nodes, four exterior LINE_2 elements and the same
Physical Groups: inflow, outflow, symmetry, farfield and fluid_domain.

The native C++ Case uses OneFLOW's Grid task. The Python reference uses the
independent STL-to-Gmsh implementation. The development runner stages both
outside the source tree; generated files do not belong in these directories.
