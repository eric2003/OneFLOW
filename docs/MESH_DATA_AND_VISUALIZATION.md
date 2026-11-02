# Mesh Data, Reconstruction, and Visualization

> This legacy copy is retained for repository history and discovery. The curated design pages under `doc/source/` are the authoritative documentation; please update those pages rather than maintaining parallel copies here.

- [Mesh architecture overview](../doc/source/mesh_architecture.rst)
- [Grid topology and validation](../doc/source/grid_topology.rst)
- [Concrete grid data contracts](architecture/grid-data-contracts-and-visualization.md)

The earlier versions of this document and the root-level `MESH_DATA_AND_VISUALIZATION.md` used different information-level taxonomies. These were design drafts, not a serialized interchange schema. The curated pages now use one shared vocabulary and distinguish current implementation facts from proposed adapters.

The detailed data-contract note remains useful for the implemented `GridSummary` / `oneflow.grid-summary` contract and the separate runtime `oneflow.mesh-info` output. Neither summary format is a reconstructable mesh package. Any richer geometry/topology package must define its own versioned schema and prove its guarantees with semantic round-trip tests.
