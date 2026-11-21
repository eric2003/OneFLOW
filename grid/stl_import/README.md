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


## Stage B1-B3: cut the analytic cylinder

The first cutting implementation deliberately starts from an explicit
analytic definition rather than the STL manifest or patch indices. The axis
is X, the centerline is y=z=0, and the lateral surface satisfies:

    y^2 + z^2 = R^2

A plane x=c intersects the cylinder in the circle:

    (x,y,z) = (c, R cos(theta), R sin(theta))

For the standard example, choose x=-1.5 and x=+1.5. Retain only the lateral
surface between those planes. Do not create end caps. This produces one
connected open surface with exactly two boundary loops.

On Windows PowerShell, run:

    python grid/stl_import/analytic_cylinder.py --cut-start -1.5 --cut-end 1.5 --output-stl "$env:USERPROFILE/.oneflow/grid-development/stl-cutting/cylinder_internal_flow/cylinder_cut.stl" --report "$env:USERPROFILE/.oneflow/grid-development/stl-cutting/cylinder_internal_flow/cut_report.json"

Then run the generic Stage A inspector independently:

    python grid/stl_import/inspect.py --input "$env:USERPROFILE/.oneflow/grid-development/stl-cutting/cylinder_internal_flow/cylinder_cut.stl" --output "$env:USERPROFILE/.oneflow/grid-development/stl-cutting/cylinder_internal_flow/stage_a_report.json"

On Linux, replace the output root with
$HOME/.oneflow/grid-development/stl-cutting/cylinder_internal_flow/.

The default discretization is 48 circumferential by 12 axial divisions. The
cut report records the analytic parameters. The separate Stage A report
re-reads only the generated STL and verifies two boundary loops, 48 edges per
loop, one connected surface component, and no non-manifold edges. The open
surface is intentionally not volume-mesh-ready.

This is a bounded analytic reference implementation, not yet a generic STL
clipper. Generic plane/triangle clipping on an unknown curved pipe is deferred
until this analytic case is proven.


## Analytic cylinder cut-to-CFD workflow

The Mesh Generation Workbench exposes a dedicated **Cylinder Cutting & CFD** entry for the `geometry.cylinder.internal_flow` asset. The user selects two axial planes and the surface discretization. The generated open lateral STL is inspected by the standalone Stage A inspector; it must contain two simple boundary loops and no non-manifold edges.

The volume mesh is not generated from that open STL. Instead, the same explicit analytic cylinder parameters are used to construct the truncated OpenCASCADE cylinder volume in Gmsh. Gmsh assigns `wall`, `inflow`, `outflow` and `fluid_domain` Physical Groups, and the MSH4 output is checked before it can be handed to OneFLOW. This separation is deliberate: an open lateral surface is not itself a closed fluid volume.

All generated artifacts are stored under `.oneflow/mesh-generation/cylinder_internal_flow/cut_workflow` by default. The **Run OneFLOW** action opens the existing Mesh-to-CFD Validation path for the current MSH: MSH staging, CGNS conversion, OneFLOW grid ingestion and a short solver smoke test. The currently registered CFD consumer is a smoke template only; it is not a cylinder-specific physical-accuracy validation.
