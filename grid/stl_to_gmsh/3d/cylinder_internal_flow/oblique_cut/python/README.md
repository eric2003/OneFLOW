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
python grid/stl_to_gmsh/run_cylinder_oblique_cut.py --gmsh gmsh
```

The runner reads `script/cutting.json`, calls the generic STL plane cutter,
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

This case validates geometry cutting and mesh construction. It is not yet a
physical-accuracy validation of an oblique-inlet CFD solution.
