# STL file I/O and compatibility cases

This directory owns STL file-format handling and focused compatibility cases.
It does not own planar topology reconstruction or STL-to-Gmsh conversion.

## Responsibilities

- `reader.py`: ASCII and binary STL record parsing. Facet order and
  coordinates are preserved; parsing does not weld vertices or infer topology.
- `surface.py`: geometric operations on triangle surfaces, including
  tolerance-based vertex identity and extraction of a single planar boundary.
  These operations do not depend on Gmsh.
- `inspect.py`: standalone file diagnostics (format, facet count, coordinate
  bounds and XY planarity). It does not generate meshes or write into the
  source tree unless the caller explicitly chooses that output path.
- `run_case.py`: Python reference runner for the STL reader compatibility
  Cases. It writes reports only inside the externally staged Case directory.
- `ascii_*`, `binary_*`, and other named directories: small, focused input
  Cases. Native C++ Cases exercise the production OneFLOW path; Python
  references provide separate parser/geometry evidence where registered.

## Inspection

```sh
python grid/stl_import/inspect.py --input path/to/model.stl --output path/to/report.json
```

The JSON report contains the detected representation, triangle count,
coordinate bounds and whether all vertices lie in a constant-Z XY plane.
Planar boundary reconstruction and meshing are deliberately outside this
command's scope.

## Python reference Cases

Run a catalogued Case through the staging-aware runner:

```sh
python grid/stl_import/run_case.py --case-directory <staged-case-directory>
```

Do not pass a checked-in source Case directory directly: generated results
must remain outside the repository.

## Test data ownership

Runtime and development Cases keep their own `input/model.stl` so each Case
remains independently runnable and understandable. Test-only compatibility
data belongs under `grid/tests/fixtures/`, not in this production STL I/O module.
The legacy square boundary-sidecar fixture used by the STL-to-Gmsh converter
tests is maintained at
`grid/tests/fixtures/stl_to_gmsh/legacy_square/`.

## Related conversion workflow

The complete STL-to-Gmsh workflow, including conversion code, sidecar layouts,
Gmsh execution and MSH validation, lives in `grid/stl_to_gmsh/`. Its Python
implementation reuses the STL reader and surface operations from this
package. See that directory's README for the supported Case progression.
