# Minimal Cube Mesh Case

This is a runnable OneFLOW case template, organized like
`grid/minimal_rectangle`. It includes the runtime control file and all
parameter files required by `OneFLOW.exe r <case-directory>`.

## Case files

- `script/control.txt` loads the five parameter files.
- `script/simu.txt` selects the grid task.
- `script/dimension.txt` selects 3D.
- `script/grid.txt` selects classic, layout-driven generation and configures
  the generated-grid, boundary, OneFLOW-grid, and summary output paths.
- `script/part.txt` configures the single-process partition defaults.
- `grid/minimalCubeLayout.txt` describes one unit cube.

The cube has eight points, twelve edges, six faces, and one block. Points 1-4
are the lower square and points 5-8 are the upper square. Lines 1-8 are the
square perimeters; lines 9-12 connect corresponding corners. Each `Add L2F`
group assigns four ordered edges to a face, and `Add F2B` assigns the six
faces to the six local block-face positions. Each edge has three nodes, so the
intended structured result is 3 x 3 x 3 nodes and 8 hexahedral cells.

## Run on Windows

Copy the complete `grid/minimal_cube` directory to a writable workspace,
for example `D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_cube`.
Keep the repository checkout separately for the executable and shared runtime
resources.

In Visual Studio, use these command arguments:

```text
r D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_cube\\
```

Or run from PowerShell:

```powershell
.\\build\\bin\\Debug\\OneFLOW.exe r D:\\OneFLOWWorkspace\\workdir2026\\grid\\minimal_cube\\
```

Run with the OneFLOW executable's directory as the working directory when
shared `system` resources are resolved relative to the executable. The case
directory must contain `script/control.txt`, all files listed by that control
file, `grid/minimalCubeLayout.txt`, and the `results` directory.

## Expected first result

The generated structured grid should contain one zone with dimensions
`ni = nj = nk = 3`, 27 nodes, coordinates spanning `[0, 1]` on each axis,
and 8 hexahedral cells. If generation fails, use the first reported runtime
error to inspect the case configuration and face ordering before changing
general mesh-generation code.

The existing mesh package exporter is 2D-only, so the rectangle's SVG/VTU
export workflow is not enabled for this 3D example yet.
