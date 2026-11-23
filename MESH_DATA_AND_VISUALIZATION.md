# Mesh Data, Reconstruction, and Visualization

> This file is retained as a legacy entry point. The curated mesh architecture documentation lives under `doc/source/`; do not add new authoritative design material here.

- [Mesh architecture overview](doc/source/mesh_architecture.rst): information levels, canonical versus derived data, snapshots, adapters, and validation strategy.
- [Grid topology and validation](doc/source/grid_topology.rst): dimensional topology, boundary ghost-cell normalization, adjacency rebuild rules, and topology invariants.
- [Concrete grid data contracts](docs/architecture/grid-data-contracts-and-visualization.md): the implemented L0 `GridSummary` contract, `oneflow.mesh-info`, and the current JSON adapter boundary.

The earlier versions of this note used a different information-level numbering from the curated architecture pages. Those labels were design drafts, not a serialized schema. Use the versioned schema identifiers in implementation and tests; do not infer a package contract from the historical L0-L4 labels.

New architecture changes should be made in the curated pages and their linked engineering records. Keep the mesh core independent of serialization and visualization libraries, and distinguish summary metadata from reconstructable geometry and topology.
