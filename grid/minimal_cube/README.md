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

## Layout relation contract

The layout file describes topology before it describes a mesh. Entity IDs in
the text file are one-based; internal vector indices are zero-based.

- `Point(id)` defines a geometric vertex.
- `Line(id)` references its two endpoint point IDs.
- `Add L2F faceId position lineId` means "add this line to this face at
  this boundary position". The four lines of a face must be ordered around
  the face boundary, with adjacent lines sharing endpoints. The position is
  one-based and defines the boundary traversal order; it is not a line ID.
- `Add F2B blockId position faceId` means "attach this face to this block's
  local face position". For a 3D block, positions 1 through 6 map to the six
  local sides. The position is one-based in the file and is converted to a
  zero-based index internally. It is not the face ID.
- `Dim`, `Ds`, and `Boundary` provide edge discretization and boundary
  metadata; they do not replace the topology relations above.

For this single-block cube, block 1 references faces 1-6 exactly once. The
face-to-block association and the local face position are separate concepts:
the former identifies which face entity is attached, while the latter tells
the block which side it occupies.

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
