# Grid Conversion Cases

This directory groups runnable conversion examples by the mesh being
demonstrated, not by the converter used to process it.

## Organization

- cases/<case-id>/ contains a self-contained, independently runnable mesh
  Case: its source mesh, OneFLOW parameter files, expected topology, and
  instructions.
- ../interoperability/conversion_matrix.json records supported format routes
  and the runnable Cases that provide concrete coverage.
- ../case_catalog.json registers each Case for GUI and command-line discovery.
- doc/source/grid_conversion_cases.rst is the documentation entry point.

A case ID describes geometry or topology (for example, two_quads or mixed_2d).
It must not include a route prefix such as gmsh_to_cgns_. The same mesh may
eventually be exercised by multiple conversion routes without copying or
renaming the geometric example.

Each Case owns its input mesh and runtime configuration. The current OneFLOW
parameter file selects a source/target pair through sourceGridType and
targetGridType. Route selection is explicit and reproducible today; a future
GUI route selector can stage a route profile into the isolated run workspace
without changing the source Case.

Conversion Cases are manually launched from the GUI Case Workbench. They are
not automatically executed by CTest, Python unittest, or GitHub Actions.
Generated files belong in the external run workspace or a Case-local results
directory when run manually; do not commit generated artifacts.
