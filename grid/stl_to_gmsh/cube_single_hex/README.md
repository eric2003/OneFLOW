# STL to Gmsh: 3D Cube with One HEXA_8 Cell

## Purpose and scope

This is the smallest volume-mesh prototype in the STL-to-Gmsh family. It
demonstrates how a closed triangulated surface is interpreted as a six-face
solid, how Gmsh constructs a volume from that surface, and how six semantic
boundary patches survive the MSH 4.1 → CGNS → OneFLOW conversion route.

The supported input is deliberately restricted to one axis-aligned cube:
12 STL triangles, eight unique corners, two coplanar triangles per square
face, equal edge lengths, and a closed two-manifold surface. This is not a
general STL-to-volume mesher.

## Input and geometric reconstruction

The C++ and Python Cases use the same cube STL input. The construction process
validates the input before invoking Gmsh:

1. weld coincident STL vertices and recover the eight bounding-box corners;
2. verify the bounding box is axis-aligned and has equal X/Y/Z edge lengths;
3. verify the 12 facets form six planar square faces with two triangles per
   face;
4. verify each undirected surface edge has exactly two incident triangles;
5. reconstruct six planar Gmsh surface entities and close them into one
   volume.

These checks are intentionally stronger than merely checking that the STL
file contains 12 facets. A malformed or differently triangulated object must
not be silently interpreted as this cube fixture.

## Mesh layout and boundary semantics

`script/layout.json` declares `mesh_dimension = 3`, the `fluid_domain`
volume Physical Group, and the semantic face-to-boundary mapping:

| Geometric face | Physical boundary | Intended meaning |
|---|---|---|
| XMin | `inflow` | inlet patch |
| XMax | `outflow` | outlet patch |
| YMin | `wall` | wall patch |
| YMax | `farfield` | far-field patch |
| ZMin | `symmetry` | symmetry plane |
| ZMax | `symmetry` | symmetry plane |

The two Z faces share a CGNS boundary-condition type but remain distinct
Physical Groups and retain independent face lists. The mapping records
semantic identity; it does not supply inlet/outlet state values or prove a
solver boundary treatment.

## Expected topology

The transfinite/recombined volume mesh has this exact acceptance contract:

| Quantity | Expected |
|---|---:|
| Nodes | 8 |
| Boundary QUAD_4 elements | 6 |
| Volume HEXA_8 cells | 1 |
| Surface Physical Groups | 6 |
| Volume Physical Groups | 1 |

The Python reference checks the same geometric and topological invariants.
C++ and Python need not assign identical node or element IDs, but must preserve
the same semantic boundary names and expected element families.

## Conversion and runtime workflow

Both implementations are registered in `grid/case_catalog.json` as 3D
STL-to-Gmsh Cases. They appear under 3D → STL to Gmsh in Run Development
Cases, Python reference cases, and Mesh-to-CFD Validation.

The validation graph checks:

1. STL → Gmsh generation and MSH 4.1 topology;
2. Gmsh → CGNS conversion;
3. CGNS → OneFLOW grid import, boundary-region map, and grid metadata;
4. a short dimension-adjusted OneFLOW solver smoke test.

For the final smoke test, the existing consumer parameter template is reused
and its runtime dimension is overridden to 3. This verifies that the imported
grid can enter the runtime path. It does not make the one-cell cube a useful
flow-resolution mesh, does not establish a physical cube-flow setup, and is
not a numerical-accuracy or convergence test.

## Tests and reproducibility

The producer tests live in the repository test hierarchy, not in this Case
directory. The Mesh-to-CFD Validation Workbench stages the Case and all
generated files under the user's external `.oneflow/mesh-validation`
workspace. Do not run a checked-in Case in-place.

The family-level runbook in `../README.md` documents the development runner,
shared artifact handoff, supported limitations, and test commands. The
general [Development Case Methodology](../../../doc/source/development_case_methodology.rst) describes the reusable process
for designing and verifying 2D and 3D Cases.

## Limitations

The current implementation does not claim support for arbitrary closed STL
surfaces, non-axis-aligned cubes, unequal edge lengths, cavities, multiple
volumes, curved CAD surfaces, or general tetrahedral/hexahedral volume
meshing. Each of those capabilities requires a separate Case and a separate
acceptance contract.
