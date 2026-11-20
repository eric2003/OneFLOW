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

## Python reference tests

```sh
python -m unittest grid.tests.python.stl_import.test_surface_validation
```

Do not run a checked-in Case in place; generated results must remain outside
the repository.
