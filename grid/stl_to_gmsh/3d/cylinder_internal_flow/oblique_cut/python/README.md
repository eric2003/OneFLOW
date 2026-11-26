# Cylinder Internal Flow — Oblique Inlet Cutting

This is a generic STL cutting reference case. The input is the same closed
cylinder STL used by the regular internal-flow cylinder case, stored as a
content-identical Git blob. It deliberately has no `.geometry.json` sidecar:
the cutter must treat `input/model.stl` as an unknown triangle shell.

## Geometry

- Cylinder axis: X
- Reference radius: 1
- Inlet plane point: `(-1.2, 0, 0)`
- Inlet unit normal: `(0.8660254037844386, 0.5, 0)`
- Inlet inclination: 30 degrees relative to the normal cross-section
- Outlet plane: `x = 1.5`, normal `(1, 0, 0)`

The inlet plane is
`cos(30°) (x + 1.2) + sin(30°) y = 0`.

Its intersection with the cylindrical wall is an ellipse. The expected
semi-axes are (R/\cos(30°)=1.1547005) and (R=1). The outlet remains a
circle of radius 1. These analytic values are verification references only;
the generic cutter must not read or use them.

## Run

From the repository root:

```console
python grid/stl_to_gmsh/3d/cylinder_internal_flow/oblique_cut/python/run_cylinder_oblique_cut.py --gmsh gmsh
```

The runner reads the shared `../script/cutting.json`, calls the generic STL plane cutter,
then sends the resulting closed shell and checksum-bound region manifest to
the existing Gmsh volume-mesh generator. Outputs are written outside the
source tree by default:

`%USERPROFILE%/.oneflow/mesh-generation/cylinder_oblique_cut/`

The generated artifacts include `wall.stl`, `inlet.stl`, `outlet.stl`,
`closed_surface.stl`, `region_manifest.json`, `cut_report.json`,
`results/cylinder_oblique_cut.msh`, and the Gmsh run report.

## Acceptance

1. The inlet boundary is measurably elliptical; its projected major/minor
   semi-axis ratio is approximately 1.1547.
2. The outlet boundary remains circular.
3. The capped shell is closed, genus zero, consistently oriented and has
   positive signed volume.
4. Gmsh produces tetrahedral volume cells and triangular boundary cells.
5. `wall`, `inflow`, `outflow` and `fluid_domain` are all non-empty
   Physical Groups, and every exterior triangle belongs to exactly one
   boundary group.

## GUI: Mesh to CFD

Open **Mesh Generation... → Oblique-cut Cylinder STL to Gmsh...**. After
**Run Cutting + Gmsh** completes successfully, choose **Convert + Run
OneFLOW...**. The GUI passes the generated
`results/cylinder_oblique_cut.msh` to the shared Mesh-to-CFD Validation
workflow; it does not regenerate the mesh. The workflow performs:

1. Existing MSH4 inspection and required Physical Group checks.
2. Gmsh MSH → CGNS conversion.
3. CGNS → native OneFLOW unstructured-grid conversion, including boundary
   region map and grid metadata.
4. A short 3-D OneFLOW solver smoke test using the registered consumer
   template.

Validation files and the run report are stored under the mesh project's
`validation/` directory. The default CFD template is a runtime smoke-test
consumer, not a physically calibrated oblique-cylinder internal-flow model.
In particular, a successful run does not establish pressure-drop, velocity
profile, or convergence accuracy for this geometry.

This case validates geometry cutting, mesh construction, format conversion,
and solver ingestion. Physical-accuracy validation of an oblique-inlet CFD
solution remains a separate milestone requiring a suitable 3-D internal-flow
CFD template and explicit inlet/outlet/wall physical data.
