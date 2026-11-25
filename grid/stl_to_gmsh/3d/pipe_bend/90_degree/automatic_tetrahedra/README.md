# 90-degree pipe bend: automatic tetrahedral mesh validation

This directory contains the expected-artifact contract and documentation for
the generated-MSH Mesh-to-CFD Validation profile. The analytic geometry
definition and its Python/C++ generators are maintained separately in the
sibling `geometry/` directory.

## GUI workflow

1. Open **Mesh Generation -> Geometry Library**.
2. Select **90° Pipe Bend — Internal Flow**, adjust the pipe/bend dimensions
   and surface discretization if needed, and generate the STL.
3. Open the generated STL in **STL to Gmsh**. The mesh project imports the
   checksum-verified geometry and region manifests together with the
   `wall`, `inflow` and `outflow` patch surfaces.
4. Generate the Gmsh volume mesh. The generic region-aware Gmsh adapter builds
   the volume and assigns dimension-2 Physical Groups `wall`, `inflow`,
   `outflow` and dimension-3 `fluid_domain`.
5. Start **Mesh-to-CFD Validation...** from the mesh project, or select
   **3D 90-Degree Pipe Bend (Generated MSH)** in the standalone 3-D validation
   tree.

The generated STL/MSH and all reports remain in the user's
`%USERPROFILE%/.oneflow/mesh-generation/` and
`%USERPROFILE%/.oneflow/mesh-validation/` workspaces. No generated geometry
or mesh fixture is stored in this source directory.

## Validation route

The Case Catalog registers this as an existing-MSH consumer profile, not as
a second mesh generator. The resolver selects the newest completed user
project whose `project.json` declares `geometry: pipe_bend_90_degree`.
It accepts only project-relative MSH paths that remain inside the
mesh-generation workspace. If no completed project is available, the GUI
reports that the generated mesh must be created first.

The validation workflow copies the MSH into its stage workspace; it does not
rerun Gmsh. It performs MSH 4 inspection, Gmsh-to-CGNS conversion,
CGNS-to-OneFLOW ingestion and the configured short CFD smoke test. Required
Physical Groups are checked before conversion:

- Dimension 2: `wall`, `inflow`, `outflow`
- Dimension 3: `fluid_domain`

For an explicitly supplied Stage Root (including CI), reports and conversion
artifacts are written directly beneath that root:

```text
oneflow-generic-curved-shell-validation/
├── mesh-validation-report.json
├── gmsh_to_cgns/
└── cgns_to_oneflow/
```

The GUI's normal project validation workspace may retain isolated run
directories so its Run History can reconstruct previous executions. That
history behavior is not applied to explicit CI/stage output directories.

The expected artifact contract is in `expected.json`. The smoke test verifies
runtime handoff only; it is not evidence of physical accuracy, pressure-loss
convergence or Dean-vortex resolution. Future physical validation should
measure pressure loss and secondary-flow structure and include mesh
refinement/convergence evidence.
