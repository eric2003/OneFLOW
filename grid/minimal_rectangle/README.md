# Minimal Rectangle Mesh Case

This case uses the standard OneFLOW runtime entry point to generate a small
layout-driven grid, convert it to OneFLOW format, write a grid-information
summary, and export SVG/VTU visualization files for inspection.

The layout describes a structured 2D rectangle, but the current OneFLOW runtime
uses unstructured grid data structures for computation. The case therefore sets
`topoType = "u"`. The structured-grid interface is reserved for future work
and is not currently the computational representation.

## Case configuration

- `script/control.txt` selects the parameter files.
- `script/simu.txt` selects `simutask = "Grid"`.
- `script/dimension.txt` selects 2D.
- `script/grid.txt` selects `gridObj = 0` (classic generation),
  `igene = 4` (layout-driven grid generation), and unstructured topology
  (`topoType = "u"`), then configures conversion from the generated grid
  to OneFLOW format.
- `grid/minimalRectangleLayout.txt` defines a unit square with three nodes
  along each edge.

The case parameters are read by OneFLOW's normal control-file loader into the
database. The generation and conversion code then reads the database through
the existing grid configuration path.

## Run and inspect on Windows or Linux

Copy this case directory to a writable workspace outside the repository. The
workspace must contain `script/control.txt`, the parameter files, and
`grid/minimalRectangleLayout.txt`. Keep the repository checkout separately;
the runner reads `case.json` and the exporter from that checkout.

From the repository root, run the same Python command on either platform:

```sh
python tools/mesh_workflow.py run \
    --case grid/minimal_rectangle/case.json \
    --workdir /path/to/workspace/grid/minimal_rectangle \
    --oneflow /path/to/build/bin/OneFLOW
```

On Windows, for example:

```powershell
python tools/mesh_workflow.py run `
    --case grid/minimal_rectangle/case.json `
    --workdir D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_rectangle `
    --oneflow build\\bin\\Debug\\OneFLOW.exe
```

Use the actual executable path produced by the local build. The runner launches
OneFLOW from the executable's directory so shared runtime resources can be
found, then checks the runtime summary, exports the generated structured grid,
and validates the exported summary and expected package files.

The PowerShell script remains as a compatibility wrapper for existing Windows
commands:

```powershell
.\\grid\\minimal_rectangle\\run.ps1 -OneFlowExe .\\build\\bin\\Debug\\OneFLOW.exe -CaseDir D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_rectangle -RepositoryRoot D:\\github\\OneFLOW
```

## Outputs and checks

All generated artifacts stay in the external case workspace:

- `results/minimalRectangle.generated.grd`: generated grid input for conversion.
- `results/minimalRectangle.inp`: generated boundary input.
- `results/minimalRectangle.ofl`: converted OneFLOW grid.
- `results/minimalRectangle.info.json`: runtime grid summary.
- `results/minimalRectangle.mesh/`: exported mesh package.

The package contains `summary.json`, `geometry.json`, `topology.json`,
`boundary.json`, `manifest.json`, `views/mesh.svg`, and `views/mesh.vtu`.
The expected counts are 9 nodes, 4 quadrilateral cells, 12 unique edges,
8 boundary edges, and 4 internal edges. The bounding box is `[0,0,0]` to
`[1,1,0]`.

Open the SVG in a browser or the VTU in ParaView. The current exporter
reconstructs 2D topology from the structured binary dimensions; it does not
yet preserve physical boundary names or general unstructured-grid semantics.
