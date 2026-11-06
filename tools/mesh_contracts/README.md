# MeshSnapshot contract

This package defines and validates the first version of the format-neutral
**interchange manifest** for mesh data. It is a tooling-side contract, not the
production C++ mesh container and not a replacement for OneFLOW's runtime
topology.

The initial contract covers:

- positive, stable node and cell IDs (source IDs need not be contiguous);
- coordinates with an explicit spatial dimension and finite numeric values;
- ordered cell-to-node connectivity with a declared element type;
- named boundary, region and interface groups with explicit dimensions;
- source-format provenance.

Connectivity order is preserved as supplied. The validator does not reorder
nodes, infer orientation, calculate face adjacency, verify geometric quality,
or prove physical boundary semantics. Group members currently refer to cell
IDs; explicit lower-dimensional boundary-entity membership is a planned
schema extension and must not be represented as if it were already supported.

Run the contract tests from the repository root:

```bash
python -m unittest tools.mesh_contracts.test_snapshot -v
```

The next production step is to implement a C++ value snapshot from the runtime
mesh and adapt Gmsh / Plot3D / SU2 readers to this contract without exposing
format-specific types in Mesh Core.
