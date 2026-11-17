# STL to Gmsh: Cube with Automatic TETRA_4 Refinement

This Case uses the same closed cube STL and semantic boundary contract as
`single_hex` and `six_tetrahedra`, but removes transfinite constraints.
Gmsh controls both surface triangulation and the tetrahedral volume mesh
using the configured characteristic `mesh_size`.

The exact node count and tetrahedron count are intentionally not prescribed.
The acceptance contract requires a non-empty TETRA_4 volume mesh, triangular
boundary facets, all six named boundary Physical Groups and the
`fluid_domain` volume group. Quadrilateral and hexahedral elements are
rejected so the Case cannot silently become a structured hex mesh.

The native C++ and Python reference Cases use the same input STL and
equivalent layout settings. Their triangulations and counts may differ across
Gmsh versions; validation therefore checks topology families and minimum
non-empty counts rather than exact connectivity. For tetrahedral output, all
exterior faces must match the emitted TRI_3 boundary elements, independent
of the number or connectivity of generated tetrahedra. See
[the cube progression](../README.md) and
[Development Case Methodology](../../../../../doc/source/development_case_methodology.rst).
