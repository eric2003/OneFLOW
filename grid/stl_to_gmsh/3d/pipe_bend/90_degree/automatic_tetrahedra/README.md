# 90-degree pipe bend: automatic tetrahedral mesh validation

This directory registers the **existing-MSH** profile for the parameterized
pipe-bend geometry. It is deliberately separate from the analytic geometry
definition in the sibling `geometry/` directory.

The actual STL and Gmsh MSH are user-owned project artifacts, not checked-in
fixtures. Generate the pipe bend from **Mesh Generation -> Geometry Library**,
open the generated STL in **STL to Gmsh**, and create the volume mesh. The
project is stored below the user's `.oneflow/mesh-generation/` workspace.

The Case Catalog exposes this profile in the 3-D **Mesh-to-CFD Validation**
tree as **3D 90-Degree Pipe Bend (Generated MSH)**. The resolver locates the
most recently modified completed project whose `project.json` declares
`geometry: pipe_bend_90_degree`; it accepts only project-relative MSH paths
that remain inside the mesh-generation workspace. If no completed project is
available, the GUI reports that the generated mesh must be created first.

The validation workflow copies the existing MSH into an isolated run
directory; it does not rerun Gmsh. A checksum-verified region manifest in the
project makes these Physical Groups mandatory before conversion:

- Dimension 2: `wall`, `inflow`, `outflow`
- Dimension 3: `fluid_domain`

The same route performs MSH 4 inspection, Gmsh-to-CGNS conversion,
CGNS-to-OneFLOW ingestion and the configured short CFD smoke test. The
smoke test verifies runtime handoff only; it is not evidence of physical
accuracy, pressure-loss convergence or Dean-vortex resolution.

The expected artifact contract is in `expected.json`. Mesh size and all
generated artifacts remain in the user's project and validation directories,
never in this source directory.
