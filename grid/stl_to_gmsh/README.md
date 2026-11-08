# Planar STL to Gmsh 2D case

This is a standalone Grid preprocessing case. It is intentionally not registered as a OneFLOW solver case and does not add any GUI functionality.

The case reads `input/body.stl`, applies the explicit physical-boundary configuration in `config/external_flow.json`, and writes a Gmsh `.geo` file and JSON report under `results/`. If a Gmsh executable is supplied, it also writes an MSH 2.2 file.

## Run

From the repository root:

```sh
python grid/stl_to_gmsh/script/run.py
```

To generate the MSH file as well:

```sh
python grid/stl_to_gmsh/script/run.py --gmsh gmsh
```

On Windows, use the same Python entry point from PowerShell and provide the installed executable name or full path with `--gmsh`.

## Boundary semantics

The STL body is an internal `Wall` boundary. The rectangular fluid domain has explicit `Inlet`, `Outlet`, `Symmetry`, and `Farfield` Physical Curve groups; the fluid region is `Fluid`. These are Gmsh physical labels, not an automatic OneFLOW solver boundary-condition mapping.

## Outputs

- `results/body.geo`
- `results/body.msh` (only with `--gmsh`)
- `results/mesh_report.json`

Generated files are ignored by Git. The checked-in case input and configuration are never modified.
