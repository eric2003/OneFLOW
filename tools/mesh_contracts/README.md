# MeshSnapshot contract

This package defines and validates the first version of the format-neutral
**interchange manifest** for mesh data. It is a tooling-side contract, not the
production C++ mesh container and not a replacement for OneFLOW's runtime
topology.

The initial contract covers:

- positive, stable node and cell IDs (source IDs need not be contiguous);
- coordinates with an explicit spatial dimension and finite numeric values;
- ordered cell-to-node connectivity with a declared element type and matching topological dimension;
- ordered face-to-node connectivity with explicit owner and optional neighbour cell IDs;
- named region groups referencing cells, and boundary/interface groups referencing faces, with explicit dimensions;
- source-format provenance.

Connectivity order is preserved as supplied. The validator does not reorder
nodes, infer orientation, calculate face adjacency, verify geometric quality,
or prove physical boundary semantics. Face records preserve the supplied node order and owner/neighbour relationship.
The validator does not yet prove that face orientation is consistent with the
solver convention, that every cell face is present exactly once, or that a
boundary face is geometrically closed. It also does not map Gmsh Physical Groups
to OneFLOW boundary conditions; it only validates the neutral group membership
contract.

Run the contract tests from the repository root:

```bash
python -m unittest tools.mesh_contracts.test_snapshot -v
```

The next production step is to implement a C++ value snapshot from the runtime
mesh and adapt Gmsh / Plot3D / SU2 readers to this contract without exposing
format-specific types in Mesh Core.
