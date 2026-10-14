# Minimal Rectangle Mesh Case

This case uses the standard OneFLOW runtime entry point to generate a small
structured grid, convert it to OneFLOW format, write a grid-information
summary, and export SVG/VTU visualization files for inspection.

## Case configuration

- `script/control.txt` selects the parameter files.
- `script/simu.txt` selects `simutask = "Grid"`.
- `script/dimension.txt` selects 2D.
- `script/grid.txt` selects `gridObj = 0` (classic generation) and
  `igene = 4` (layout-driven grid generation), then configures conversion
  from the generated structured grid to OneFLOW format.
- `grid/minimalRectangleLayout.txt` defines a unit square with three nodes
  along each edge.

The case parameters are read by OneFLOW's normal control-file loader into the
database. The generation and conversion code then reads the database through
the existing grid configuration path. The unit test does not implement this
workflow.

## Run directly

Copy this case directory to a workspace outside the repository, for example
`D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_rectangle`. Ensure the copied
directory contains `script/control.txt`, the parameter files, and the
`grid/minimalRectangleLayout.txt` input.

Run the production executable from its build directory so it can locate the
shared `system` resources:

```powershell
cd D:\\github\\OneFLOW\\build\\bin\\Debug
.\\OneFLOW.exe r D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_rectangle\\
```

The configured outputs are written under the case directory:

- `results/minimalRectangle.generated.grd`: generated structured-grid input
  for the conversion step.
- `results/minimalRectangle.inp`: generated boundary input.
- `results/minimalRectangle.ofl`: converted OneFLOW grid.
- `results/minimalRectangle.info.json`: grid metadata emitted by the
  production grid conversion path.

## Export and check visualization

After the OneFLOW run succeeds, export the generated structured grid and check
its expected geometry and topology. Run this from a checkout of OneFLOW; the
case directory can remain external:

```powershell
.\\grid\\minimal_rectangle\\run.ps1 -OneFlowExe .\\build\\bin\\Debug\\OneFLOW.exe -CaseDir D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_rectangle -RepositoryRoot D:\\github\\OneFLOW
```

The helper invokes the same `OneFLOW.exe r <case-directory>` command, then
writes `results/minimalRectangle.mesh/summary.json`, `geometry.json`,
`topology.json`, `boundary.json`, `manifest.json`, `views/mesh.svg`, and
`views/mesh.vtu`. The expected counts are 9 nodes, 4 quadrilateral cells,
12 unique edges, 8 boundary edges, and 4 internal edges; the bounding box is
`[0,0,0]` to `[1,1,0]`.

Open the SVG in a browser or the VTU in ParaView. All generated artifacts stay
in the external case workspace, not in the checked-in case template.
