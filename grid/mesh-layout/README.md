# Mesh Layout Cases

Mesh layout cases are organized by spatial dimension:

- `2d/` contains planar layouts and flat-plate templates.
- `3d/` contains volumetric layouts.

Case names describe the specific geometry or topology being demonstrated.
For two-block cases, the name identifies the shared entity and the interface
discretization. For example, `two_blocks_shared_edge_endpoints_only` is the
two-cell baseline; `two_blocks_shared_edge_one_interior_node` adds exactly
one node inside the shared edge.

Keep each example focused on one new condition. Do not combine unrelated
boundary-condition, mesh-quality, conversion, or solver demonstrations into
a foundational layout case. Add a new case when a new topology or discretization
condition needs its own test evidence.
