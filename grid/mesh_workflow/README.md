# Mesh Workflow Case

This is a standalone grid-processing case template, organized like the existing
`test/grid/rae2822convert` runtime case but kept outside the test tree.

The repository copy is a template. Run it from a separate workspace directory
so logs and generated files do not pollute the source checkout.

## Current workflow

1. Read the 2D RAE 2822 Plot3D grid and boundary-condition input.
2. Run the regular OneFLOW grid-conversion entry point (`gridObj = 1`).
3. Write the OneFLOW grid file to `results/rae2822_2d.ofl`.
4. Write a JSON grid-information summary to `results/rae2822_2d.info.json` with dimension, node/face/cell counts, boundary-face count, and coordinate bounding box. The summary contains metadata only, not mesh connectivity or coordinates.

The case configuration is in `script/grid.txt`. Input and output paths are
relative to the active case directory.

## Run from Visual Studio or a terminal

Use the actual case directory in your workspace, not this repository template.
The equivalent command-line invocation is:

```text
OneFLOW.exe r D:\OneFLOWWorkspace\workdir2026\mesh_workflow\
```

In Visual Studio 2022, set the executable's command arguments to
`r D:\OneFLOWWorkspace\workdir2026\mesh_workflow\`. Keep the working
directory set to the OneFLOW executable's installation/build directory when
the runtime needs to resolve shared `system` resources. Case-relative input,
logs, and generated files should remain under the external case directory.

The optional `run.ps1` launcher is only a convenience; it is not required.

## Preflight validation

Before launching OneFLOW, use the standard-library Python validator to check
the control-file includes, required grid-conversion parameters, source grid
and boundary input files, and output directory:

```powershell
python .\grid\mesh_workflow\tools\validate_case.py D:\OneFLOWWorkspace\workdir2026\mesh_workflow
```

Add `--create-output-dirs` if the configured output directory does not yet
exist. The validator currently supports this case's `simutask = Grid`,
`gridObj = 1` conversion workflow. It also checks that `gridLayoutFileName` is
present because `GridConfig::FromDataBase()` reads that parameter even for a
conversion; it does not require the referenced layout file to exist unless the
selected workflow actually opens it. If `sourceGridCaseDir` is set, the
validator checks source grid and boundary paths relative to that source case,
matching the grid reader's case-root selection. Source and target format names
are also checked against OneFLOW's registered format tokens and the currently
implemented source-to-target conversion pairs: Plot3D to OneFLOW/CGNS, SU2 to
OneFLOW, and CGNS to OneFLOW/CGNS. This is a preflight check, not a replacement
for OneFLOW's own parameter parser or runtime validation.

## Parameter configuration and future tools

Keep `control.txt` as the authoritative list of parameter files. Some legacy
parameters are accepted but may not be used by this particular workflow; do
not remove them until their read paths and compatibility requirements have
been traced.

The next step is to classify parameters against their actual consumers, then
introduce a small case-profile model that can generate the existing text
parameter files and drive validation. A GUI can later edit the same profile;
it should not introduce a second, independent parameter system.

## Next increments

Extend this as an end-to-end grid case to emit a mesh summary, export topology
and boundary metadata through supported grid interfaces, and generate
visualization files. Prefer exporting a snapshot from live grid objects. If a
postprocessor reads the OneFLOW grid file, it must follow the actual file layout
rather than assume a simplified binary header.
