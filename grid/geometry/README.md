# OneFLOW Geometry Assets

This directory defines the format-independent geometry asset registry. It is not an STL directory and is not a runnable Case catalog.

- `asset.schema.json` defines the metadata contract for one geometry asset.
- `catalog.json` registers the initial analytic/parametric reference assets.

## Asset versus representation

A geometry asset is the identity, design intent and provenance of a geometric object. A representation is one concrete way to express or store it. One asset may have multiple representations, for example:

- analytic equations and parameter values;
- a parametric or procedural source definition;
- exact CAD boundary representation (B-rep), such as STEP or native CAD;
- tessellated surface data, such as STL, OBJ, PLY or 3MF;
- volume meshes, such as Gmsh MSH, CGNS, VTU or solver-specific formats;
- point clouds or implicit fields.

The catalog must not imply that every format is currently readable, writable, or supported by every OneFLOW workflow. Format capability belongs to the specific importer/exporter or application adapter.

## Provenance and derivation

Every representation should be traceable to its asset definition and, where applicable, to a generator, tutorial/example, source revision and parent representation. Generated artifacts should record their parameters, units, generator version and content checksum in a sidecar manifest or equivalent metadata. A derived mesh must retain a link to the exact geometry representation from which it was created.

Historical fixtures with unknown original authorship or generator must remain marked as unknown; a newer reproducible generator must not be retroactively claimed as their source.

## Boundaries

- `grid/case_catalog.json` remains the catalog of runnable development Cases.
- `grid/geometry/catalog.json` is the catalog of reusable geometric assets.
- `grid/stl_to_gmsh/` remains a mesh-generation workflow and its examples; it is a consumer/producer of representations, not the geometry registry.
- Runtime artifacts belong in the user's OneFLOW workspace, not in this source-controlled registry.

The registry is deliberately metadata-first. It does not yet promise a generic CAD kernel, universal format conversion, or a stable cross-language geometry API. Those capabilities must be added behind explicit adapters and verified per format.
