# Mesh Workflow Case

This is a standalone grid-processing case, organized like the existing
`test/grid/rae2822convert` runtime case but kept outside the test tree.

## Current workflow

1. Read the bundled 2D RAE 2822 Plot3D grid and boundary-condition input.
2. Run the regular OneFLOW grid-conversion entry point (`gridObj = 1`).
3. Write the OneFLOW grid file to `results/rae2822_2d.ofl`.

The case configuration is in `script/grid.txt`. Input and output paths are
relative to this case directory, so the case can be copied and run independently.

## Run

From the OneFLOW executable directory, where the runtime `system` configuration
is available, invoke the executable with this case directory using the same
argument convention as the existing case runner:

```powershell
.\OneFLOW.exe 0 <absolute-path-to-repository>\grid\mesh_workflow
```

Create the `results` directory before running if it does not exist.

## Next increments

Keep this as a working case rather than a unit-test fixture. Extend it in
small end-to-end steps to emit a mesh summary, export topology and boundary
metadata through supported grid interfaces, and generate visualization files.
The visualization stage must read the actual OneFLOW grid-file layout or,
preferably, consume a snapshot exported from the live grid objects; it must not
assume a simplified binary header.
