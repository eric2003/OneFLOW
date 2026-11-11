# STL to Gmsh: Cube with One HEXA_8 Cell

This is the first 3D volume Case in the cube progression. It validates the
smallest possible hexahedral volume mesh: one closed cube, eight corner nodes,
six QUAD_4 boundary faces and one HEXA_8 volume cell.

The input is a closed axis-aligned cube STL with 12 surface triangles, two
coplanar triangles per square face, and a closed two-manifold edge incidence.
The conversion reconstructs six planar CAD surfaces and one volume. Gmsh
uses transfinite curves, surfaces and volume, then recombines the result.

The six boundary Physical Groups are `xmin_inflow`, `xmax_outflow`,
`ymin_wall`, `ymax_farfield`, `zmin_symmetry` and `zmax_symmetry`;
the volume group is `fluid_domain`.

Both native C++ and Python reference implementations are registered in the
Case catalog and appear in the 3D sections of Run Development Cases, Python
Reference Cases and Mesh-to-CFD Validation. Runtime outputs are staged outside
the repository. See [the cube progression](../README.md) and
[Development Case Methodology](../../../../doc/source/development_case_methodology.rst).

This Case demonstrates only one axis-aligned cube and one HEXA_8 cell. It is
not a general STL-to-hexahedral mesher or a numerical-accuracy benchmark.
