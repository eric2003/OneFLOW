# OneFLOW mesh interoperability

This directory is the central inventory and design entry for mesh-format
interoperability. It is intentionally separate from solver Cases and from
format-specific parser implementations.

## Target organization (introduced incrementally)

- `formats/`: canonical tiny fixtures and format-specific notes. A source
  mesh is stored once, regardless of how many conversion routes consume it.
- `routes/`: route-specific Case configuration, expected outcomes and
  generated evidence. These directories reference format fixtures; they do
  not duplicate them.
- `contracts/`: shared semantic checks for topology, identity, boundaries
  and round trips.
- `conversion_matrix.json`: machine-readable inventory of routes,
  implementation status and limitations.

These subdirectories are expanded when their first real fixture, route or
contract is migrated; they are not a requirement to create empty folders.
Existing manual Cases such as `grid/gmsh_to_cgns` remain in place during this
incremental transition. Its `grid/square.msh` is the current canonical Gmsh
fixture and must not be copied into every route.

## Architecture

```text
format reader / writer adapters
             |
             v
  neutral mesh semantics
             |
             v
      OneFLOW Mesh Core
             |
             v
       solver runtime

CGNS is a preferred on-disk interchange hub, not the runtime mesh model.
```

A composed route such as Gmsh -> CGNS -> OneFLOW must write and reopen a real
CGNS file. The interchange artifact and the native OneFLOW grid output are
separate files; the latter must never overwrite the former. Tests should
inspect the constructed OneFLOW mesh model (not merely check that a file
exists). This validates serialization, indexing conventions, the CGNS reader
and the internal topology construction.

For the current 2D Gmsh fixture, boundary validation is semantic rather than
count-only: the OneFLOW boundary record must preserve the mapped BC type and
registered physical-group name, and each pair must remain attached to its
expected geometric edge (left=inflow, right=outflow, bottom=symmetry,
top=farfield). This catches a particularly dangerous class of conversion
defects where all four boundary faces exist but their solver meanings have
been permuted.

## Current routes

See `conversion_matrix.json` for the maintained route inventory. Currently
implemented paths include Gmsh -> CGNS, Gmsh -> OneFLOW (composed through
CGNS), CGNS -> OneFLOW, Plot3D -> CGNS, Plot3D -> OneFLOW and SU2 -> OneFLOW.

## Adding a route

1. Reuse or add the smallest source-format fixture under the format owner.
2. Declare exact input scope, element families, dimensions and semantic
   limitations in the route manifest.
3. Reuse existing adapters and common contracts before adding conversion
   code.
4. Add independent round-trip tests for coordinates, ordered connectivity,
   cell/face topology, group membership and boundary meaning.
5. Keep generated files and reports under temporary test directories or
   Case `results/`; never write test artifacts into source fixtures.
