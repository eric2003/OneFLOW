# Mesh Workflow Case

This is a standalone grid-processing case template, organized like the existing
`test/grid/rae2822convert` runtime case but kept outside the test tree.

The repository copy is a template. Run it from a separate workspace directory
so logs and generated files do not pollute the source checkout.

## Current workflow

1. Read the 2D RAE 2822 Plot3D grid and boundary-condition input.
2. Run the regular OneFLOW grid-conversion entry point (`gridObj = 1`).
3. Write the OneFLOW grid file to `results/rae2822_2d.ofl`.

The case configuration is in `script/grid.txt`. Input and output paths are
relative to the active case directory.

## Run on Windows

Use the actual case directory in your workspace, not this repository template:

```powershell
.\run.ps1 -OneFlowExe "D:\OneFLOWWorkspace\build\bin\Debug\OneFLOW.exe" -CaseDir "D:\OneFLOWWorkspace\workdir2026\mesh_workflow"
```

The launcher invokes the standard multi-case mode, equivalent to:

```text
OneFLOW.exe r D:\OneFLOWWorkspace\workdir2026\mesh_workflow\
```

It runs the executable from its installation directory so OneFLOW can locate
its `system` resources, while all case-relative input, logs, and output files
remain under the external case directory. The case path is explicit and
required; the script will not silently run the copy inside the source tree.

## Parameter configuration and future tools

Keep `control.txt` as the authoritative list of required parameter files.
Unused legacy parameters may remain for compatibility, but parameters required
by the selected workflow must not be omitted.

A useful next improvement is a small parameter generator/validator that reads
a case profile, emits the existing text parameter files, and checks referenced
input files and output directories before launching OneFLOW. A GUI can be added
later as a front end to the same profile model; it should not introduce a
second, independent parameter system.

## Next increments

Extend this as an end-to-end grid case to emit a mesh summary, export topology
and boundary metadata through supported grid interfaces, and generate
visualization files. Prefer exporting a snapshot from live grid objects. If a
postprocessor reads the OneFLOW grid file, it must follow the actual file layout
rather than assume a simplified binary header.
