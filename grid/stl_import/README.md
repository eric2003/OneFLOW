# STL import and surface diagnostics

This package owns STL file-format handling and non-destructive surface
inspection. It does not own CAD reconstruction, clipping, cap generation,
CFD boundary assignment, or volume meshing.

## Import and inspection layers

- `reader.py`: ASCII and binary STL parsing. Facet coordinates are preserved;
  parsing does not weld vertices or infer topology.
- `validation.py`: scale-aware diagnostics after temporary tolerance-based
  vertex welding. It reports degenerate/duplicate facets, boundary and
  non-manifold edges, orientation conflicts, connected components, closedness,
  Euler characteristic/genus where applicable, and signed volume.
- `surface.py`: topology operations used by planar conversion.
- `inspect.py`: standalone JSON report command. It never repairs the source,
  assigns inlet/outlet/wall semantics, or generates a mesh.
- `run_case.py`: Python reference runner for focused STL reader Cases.

## Inspect an arbitrary STL

```sh
python grid/stl_import/inspect.py --input path/to/model.stl --output path/to/stl_report.json
```

The default weld tolerance is 1e-9 times the bounding-box diagonal, with a
1e-12 numerical floor. Override it with `--weld-tolerance` in model units.
STL does not encode units; the report marks units as unknown and never silently
converts coordinates.

Report status:
- `valid`: no detected degenerate/duplicate facets, non-manifold edges or
  orientation conflicts; one connected component.
- `warning`: readable and inspectable, but open or multi-component. Open
  surfaces may be intentional and are not import failures.
- `invalid`: detected degenerate/duplicate facets, non-manifold edges or
  inconsistent edge orientation need attention.

`volume_mesh_ready` is stricter than readable: it requires one closed,
edge-manifold, consistently oriented component with positive signed volume.
Genus is reported, not restricted to zero. Self-intersection testing is not
performed in this stage. No repair or CFD boundary inference is performed.

## Stage A: inspect the canonical cylinder without patch metadata

The generic inspector is also used as a deliberately blind baseline for the
internal-flow cylinder. It receives only the checked-in binary STL; it does
not open the sibling `model.geometry.json`, use patch indices, infer inlet or
outlet semantics, call Gmsh, or modify the source geometry.

On Windows PowerShell:

```powershell
python grid/stl_import/inspect.py --input grid/stl_to_gmsh/3d/cylinder_internal_flow/automatic_tetrahedra/python/input/model.stl --output "$env:USERPROFILE/.oneflow/grid-development/stl-inspection/cylinder_internal_flow/stage_a_report.json"
```

On Linux:

```sh
python grid/stl_import/inspect.py --input grid/stl_to_gmsh/3d/cylinder_internal_flow/automatic_tetrahedra/python/input/model.stl --output "$HOME/.oneflow/grid-development/stl-inspection/cylinder_internal_flow/stage_a_report.json"
```

The regression test runs this same CLI against the checked-in STL and checks
that the surface is one closed, consistently oriented, genus-zero component
with no boundary or non-manifold edges. The report is independent of CFD patch
classification; Stage B will introduce user-specified cutting planes later.

## Python reference tests

```sh
python -m unittest grid.tests.python.stl_import.test_surface_validation
```

Do not run a checked-in Case in place; generated results must remain outside
the repository.
