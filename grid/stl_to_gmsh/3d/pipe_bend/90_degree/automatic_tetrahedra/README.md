# 90-degree pipe bend: generated tetrahedral mesh Case

This is a runnable 3-D STL-to-Gmsh producer Case, registered in
`grid/case_catalog.json` as `stl_to_gmsh_pipe_bend_90_degree`. In the
Mesh-to-CFD Validation dialog it appears under **3D → STL to Gmsh**, not under
Existing Gmsh Mesh.

The Case runner performs the complete mesh-production stage:

1. Generate the parameterized 90-degree circular pipe-bend STL.
2. Write the `wall`, `inflow`, and `outflow` region STL files and the
   checksum-linked `region_manifest.json` / `.geometry.json` lineage.
3. Invoke the shared region-aware Gmsh tetrahedral volume mesher to create
   `results/pipe_bend_90_degree.msh`.
4. Validate MSH4 and require dimension-2 `wall`, `inflow`, `outflow`
   and dimension-3 `fluid_domain` Physical Groups.
5. Continue through Gmsh-to-CGNS, CGNS-to-OneFLOW, and the configured short
   CFD smoke test.

The geometry generator and volume mesher remain reusable components; this
Case runner is only their workflow adapter. No generated STL or MSH is stored
in the source tree. Case outputs are written to the user's OneFLOW workspace.

The CFD smoke test verifies runtime handoff only. It is not evidence of
physical accuracy, pressure-loss convergence, or Dean-vortex resolution.
Future physical validation should add pressure-drop measurements, secondary
flow analysis, and mesh refinement/convergence evidence.
