# Two Connected Squares Mesh Case

This case is the 2D counterpart to `two_connected_cubes`. It generates two
adjacent unit squares that share one edge, using the normal OneFLOW grid task.

The squares occupy:

- Block 1: `[0, 1] x [0, 1]`
- Block 2: `[1, 2] x [0, 1]`

The shared interface is the segment `x = 1`. Each outer edge and the shared
edge are discretized with three points, so each square contains four quadrilateral
cells. The expected mesh has 15 unique nodes, 8 cells, 21 unique edges, 12
boundary edges, and 9 internal edges.

## Parameter migration

`case.json` is the versioned case profile and includes the current grid
generation options, output names, and expected geometric/topological counts.
The existing `script/*.txt` files remain the runtime-compatible inputs for
this migration stage. The values in `script/grid.txt` intentionally mirror
the profile's `gridOptions`; OneFLOW still reads the text files directly.

This is an incremental migration point, not a runtime switch to JSON. Later
work can generate the legacy text files from the profile and add consistency
checks before moving other cases. Older cases and parameter files remain valid.

The explicit IO roots use `case`, so input and output paths resolve relative to
the active case directory. Generated files remain under `results/`.

## Run

Copy this directory to a writable workspace outside the repository and run:

```text
OneFLOW.exe r <case-directory>
```

For example:

```powershell
.\build\bin\Debug\OneFLOW.exe r D:\OneFLOWWorkspace\grid\two_connected_squares\
```

Use the executable path from the local build. Run from the executable/build
directory if shared runtime resources are resolved relative to it.

## Expected checks

The generated summary should report dimension 2, two zones, 15 nodes, 8
quadrilateral cells, 21 unique edges, 12 boundary edges, and 9 internal edges.
The bounding box is `[0, 0, 0]` to `[2, 1, 0]`. The shared interface must be
present once in the combined topology and attached to both blocks.

The profile can be checked with:

```sh
python tools/mesh_cases.py validate two_connected_squares
python tools/test_mesh_cases.py
```
