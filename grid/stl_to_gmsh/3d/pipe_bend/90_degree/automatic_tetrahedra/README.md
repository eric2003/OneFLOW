# 90-degree pipe bend: incremental 3-D Case

This is the first curved 3-D pipe Case being delivered through the shared
OneFLOW Case lifecycle. It is a normal Case, not a pipe-bend-specific workflow
framework. Its lifecycle contract is
`case-lifecycle.json`; the Case Catalog entry is
`stl_to_gmsh_pipe_bend_90_degree`.

## Implementations

The same Case offers two geometry implementations:

- `python_reference`: Python reference geometry producer.
- `native_cpp_geometry`: C++17 geometry producer.

Both produce the same semantic artifacts and region identities. Both hand off
to the shared, region-aware Gmsh volume mesher. The C++ option is a native
geometry implementation; its current adapter still uses the shared Python
Gmsh service. It is not a native OneFLOW solver.

## Current producer milestone

The current verification scope is geometry, semantic surface regions, and
Gmsh volume-mesh generation:

1. Generate the parameterized 90-degree circular pipe-bend surface.
2. Export the closed fluid-domain STL, the `wall`, `inflow`, and `outflow`
   patch STL files, and checksum-linked geometry/region manifests.
3. Validate that the named patches form an exact, non-overlapping partition
   of the closed surface.
4. Invoke the shared region-aware Gmsh tetrahedral mesher.
5. Validate MSH4 structure, first-order tetrahedral volume cells, boundary
   connectivity, and the required Physical Groups:
   dimension-2 `wall`, `inflow`, `outflow`, and dimension-3 `fluid_domain`.

The producer runners write generated inputs and results to the user's OneFLOW
workspace. They do not put generated STL/MSH files in the source tree. The
existing standalone producer entry points are workflow adapters; geometry
construction and volume meshing remain reusable services.

## Remaining lifecycle stages

The MSH4 inspector calculates scale-independent TETRA_4 mean-ratio
quality, minimum and mean element quality, P05/P50/P95 quality percentiles,
counts below the 0.1/0.2/0.3 diagnostic levels, minimum and total cell volume,
and the ten lowest-quality elements (element/node tags and centroid). The
distribution and worst-element records are persisted in mesh-validation
evidence; scalar statistics are also exposed as run metrics and restored during
evidence recovery. The 0.1/0.2/0.3 levels are diagnostic buckets only, not
acceptance limits. Degenerate tetrahedra are rejected. The quality stage is
not yet an acceptance gate: no pass/fail thresholds are applied until
case-appropriate limits and their rationale are declared.

The following stages are declared in the lifecycle contract but are not yet
delivered as verified stages:

- **Mesh quality acceptance:** declare and apply justified tetrahedral quality
  thresholds, then verify the actual pipe-bend mesh against them.
- **MSH to CGNS:** use the shared converter and verify cell counts, dimensions,
  connectivity, and named boundary preservation.
- **OneFLOW ingestion:** read the CGNS grid into the native grid model and
  verify inlet, outlet, wall, dimensions, and runtime environment evidence.
- **CFD smoke:** execute a bounded solver run and verify finite residuals and
  required result artifacts.
- **CFD physical validation:** independently validate pressure loss and
  secondary-flow behavior, with mesh-refinement/convergence evidence.

A stage is not considered passed merely because its implementation exists.
The lifecycle contract keeps implementation status separate from verification
status, and the run evidence is authoritative. The short CFD smoke test, when
added, will prove runtime handoff only; it will not prove physical accuracy,
pressure-loss convergence, or Dean-vortex resolution.

## Reuse requirement

Future S-shaped pipes, other bend angles, branched ducts, and imported curved
surfaces must use the same lifecycle vocabulary and shared geometry/meshing/
conversion services. They must not require another hard-coded per-case
workflow runner. Case-specific geometry parameters and semantic region
definitions belong in the Case inputs and contract.
