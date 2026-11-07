# Gmsh minimal example matrix

This directory is reserved for minimal MSH fixtures grouped by the contract
they exercise. Fixtures should be deliberately tiny: the purpose is to verify
format structure, topology, ordered connectivity, and conversion semantics,
not solver convergence or mesh performance.

Planned fixture groups:

- `read_only/`: a single TRI_3 mesh with sparse node tags; then one fixture
  per parser feature such as Physical Groups and entity ownership.
- `convert_2d/`: the four-node square with TRI_3 cells and four named
  boundary edges; a separate QUAD_4 fixture when mixed/topological coverage
  requires it.
- `convert_3d/`: one-cell TETRA_4, HEXA_8, PENTA_6, and PYRA_5 meshes, each
  with the minimum codimension-one faces needed to verify boundary sections.
- `higher_order_read_only/`: minimal TRI_6 and TETRA_10 fixtures to verify
  node counts and ordered high-order connectivity before writer support.
- `invalid/`: malformed node references, duplicate IDs, and unsupported
  element families with explicit expected diagnostics.

The current canonical 2-D boundary conversion fixture remains
`../grid/square.msh`. Do not duplicate it here. Add a fixture only when a
focused reader or CGNS round-trip test consumes it. Keep automated tests
independent of the manual GUI Case runner: parser tests need no CGNS library,
conversion tests use temporary output files and the CGNS read API, and the
manual Case remains responsible for the full workflow and reports.
