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

## Current verified workflow scope

The lifecycle contract now covers the complete shared execution chain, from
geometry generation through a bounded native CFD smoke run:

1. Generate the parameterized 90-degree circular pipe-bend surface.
2. Export the closed fluid-domain STL, the `wall`, `inflow`, and `outflow`
   patch STL files, and checksum-linked geometry/region manifests.
3. Validate that the named patches form an exact, non-overlapping partition
   of the closed surface.
4. Invoke the shared region-aware Gmsh tetrahedral mesher and validate the
   MSH4 topology, boundary connectivity, and Physical Groups:
   dimension-2 `wall`, `inflow`, `outflow`, and dimension-3 `fluid_domain`.
5. Convert MSH4 to CGNS using the shared OneFLOW converter.
6. Read CGNS into the OneFLOW native grid and retain the boundary-region map
   and grid metadata. In particular, `wall`, `inflow`, and `outflow`
   must remain identifiable through the conversion handoff.
7. Run the shared CFD consumer with a bounded step count (the current smoke
   invocation uses two steps) and verify successful process completion and
   the expected residual/restart/visual output artifacts.

The Mesh-to-CFD Validation report is authoritative for whether each stage
actually ran and passed. A completed Gmsh-to-CGNS conversion alone does not
prove OneFLOW ingestion or CFD execution. The shared workflow writes each
stage's inputs, results, and logs under the external project `validation/`
directory; generated data is not written into the source tree.

## Mesh quality and physical validation

The MSH4 inspector reports TETRA_4 mean-ratio quality, minimum and mean
quality, P05/P50/P95 percentiles, counts below diagnostic levels, minimum and
total cell volume, and worst-element locations. These are advisory diagnostics
and do not block conversion, ingestion, or smoke execution. The current
smoke run is a software/runtime handoff check, not a physical validation.

The following remain outside the current verified scope:

- **Mesh quality acceptance:** define justified thresholds and test mesh
  sensitivity before claiming production suitability.
- **CFD physical validation:** compare pressure loss and secondary-flow
  behavior with an independent analytical, experimental, or published
  reference, including mesh-refinement/convergence evidence.

The two-step CFD smoke test does not prove pressure-loss accuracy, convergence,
or Dean-vortex resolution. Stage implementation status and verification
status remain separate; only actual run evidence can establish a passed stage.

## Reuse requirement

Future S-shaped pipes, other bend angles, branched ducts, and imported curved
surfaces must use the same lifecycle vocabulary and shared geometry/meshing/
conversion services. They must not require another hard-coded per-case
workflow runner. Case-specific geometry parameters and semantic region
definitions belong in the Case inputs and contract.
