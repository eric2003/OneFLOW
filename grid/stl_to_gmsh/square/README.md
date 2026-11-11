# Square STL mesh variants

The square geometry is shared by the single-quad and two-triangle Cases.
The STL is the same two-facet representation in both Cases; its diagonal is
removed when extracting the outer boundary and is not itself a CFD cell
interface.

| Variant | Layout | Surface recombination | Output |
|---|---|---|---|
| square/single_quad/cpp and square/single_quad/python | 1 x 1 | enabled | 1 QUAD_4 |
| two_triangles/cpp and two_triangles/python | 1 x 1 | disabled | 2 TRI_3 |
| automatic_triangles/cpp and automatic_triangles/python | automatic | not used | unstructured TRI_3 controlled by mesh_size |

The fixed-topology variants have four nodes, four exterior LINE_2 elements and the same
Physical Groups: inflow, outflow, symmetry, farfield and fluid_domain. The
`automatic_triangles` pair omits `transfinite_cells` and lets Gmsh choose the
unstructured triangle topology from `mesh_size`; the element count is not fixed.

The native C++ Case uses OneFLOW's Grid task. The Python reference uses the
independent STL-to-Gmsh implementation. The development runner stages both
outside the source tree; generated files do not belong in these directories.


## Case authoring and execution notes

The three square variants are a deliberately incremental teaching sequence:

1. `single_quad`: verify that the two-facet STL is interpreted as one planar
   surface and recombined into one QUAD_4 cell.
2. `two_triangles`: disable recombination and use the same boundary to obtain
   the two TRI_3 cells separated by the square diagonal.
3. `automatic_triangles`: remove transfinite constraints and let Gmsh choose
   an unstructured TRI_3 subdivision from `mesh_size`.

For the third Case, `mesh_size` is a target characteristic length, not an
exact element count. The C++ and Python implementations use byte-identical
STL and equivalent layout inputs, but their generated triangulations need not
have identical connectivity. See
`../automatic_triangles/README.md` for the complete input-to-evidence
walkthrough and [Development Case Methodology](../../../doc/source/development_case_methodology.rst) for the reusable authoring
standard.
