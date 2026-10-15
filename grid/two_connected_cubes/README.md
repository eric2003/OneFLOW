# Two Connected Cubes Mesh Case

This case extends `minimal_cube` to two adjacent 3D blocks. It is a readable,
modifiable input template for the standard OneFLOW entry point:

```text
OneFLOW.exe r <case-directory>
```

Copy the complete case directory to a writable workspace outside the repository
before running it. The repository directory contains inputs only; generated
files belong under the workspace's `results/` directory.

## Geometry and topology

The two unit cubes occupy:

- Block 1: `[0, 1] x [0, 1] x [0, 1]`
- Block 2: `[1, 2] x [0, 1] x [0, 1]`

They share the plane `x = 1`. The layout uses 12 unique corner points, 20
unique edges, 11 faces, and two blocks. The shared face is represented once in
the face topology and attached to both blocks. Its boundary marker is -1 in
the layout; the remaining faces are exterior faces.

Every edge has two nodes in this minimal case, so each block generates one
hexahedral cell. The expected structured source representation has two zones,
each with dimensions `2 x 2 x 2`, and four matching nodes on the shared plane.
The current runtime conversion path may store the resulting computational
grid using unstructured data structures; this case specifically exercises the
structured multi-block layout before conversion.

## Case files

- `script/control.txt` is the entry point and lists the parameter files.
- `script/simu.txt` selects the grid task.
- `script/dimension.txt` selects 3D.
- `script/grid.txt` selects layout-driven generation and names the outputs.
- `script/part.txt` supplies the standard single-process partition defaults.
- `grid/twoConnectedCubesLayout.txt` defines points, edges, faces, blocks,
  edge discretization, and boundary markers.
- `results/.gitkeep` keeps the expected output directory in the case template.

The layout relation semantics are the same as in `minimal_cube`:
`Add L2F` orders edges around a face, while `Add F2B` assigns a face to a
block-local side. Entity IDs and local side positions in the file are
one-based.

## Run

Copy `grid/two_connected_cubes` to an external workspace, then launch OneFLOW
with that copied directory as the case argument. Keep the process working
directory at the OneFLOW executable/build directory if shared runtime
resources are resolved relative to the executable.

For example, in Visual Studio 2022, set the command arguments to:

```text
r D:\OneFLOWWorkspace\grid\two_connected_cubes\
```

Or from PowerShell:

```powershell
.\build\bin\Debug\OneFLOW.exe r D:\OneFLOWWorkspace\grid\two_connected_cubes\
```

Use the actual path to your built executable.

## Expected checks

After generation, inspect `results/twoConnectedCubes.generated.grd`,
`results/twoConnectedCubes.inp`, `results/twoConnectedCubes.ofl`, and
`results/twoConnectedCubes.info.json`. The summary should report dimension
3 and the bounding box `[0, 0, 0]` to `[2, 1, 1]`.

The key geometry/topology invariants are:

- two zones / blocks;
- one hexahedral cell per block;
- four coincident shared-interface nodes, with matching coordinates;
- no gap or overlap at `x = 1`;
- exterior boundary faces remain distinct from the shared interface.

This case is intentionally small. Use failures here to improve the layout
contract and generation pipeline before generalizing to larger multi-block
meshes.
