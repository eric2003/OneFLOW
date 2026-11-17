# Shared mesh interoperability contracts

Contracts describe semantic invariants independently of a particular
conversion route. Add a contract only when it can be reused or provides a
clear independent validation boundary.

Initial contract families:

- Topology: dimensions, supported element family, node arity and ordered
  connectivity.
- Identity: sparse source IDs, entity scopes and explicit target index maps.
- Boundaries: region identity, membership, orientation and explicit semantic
  mappings.
- Round trips: coordinates, cell/face topology and metadata preserved across
  serialization and re-import.

The tooling-side MeshSnapshot v1 validator is under
`tools/mesh_contracts`. It is not a replacement for the production C++
Mesh Core.
