# STL format inspection and compatibility cases

Each STL compatibility case is organized by problem, with separate `cpp/` and
`python/` implementations when both exist.

```text
grid/stl/
  ascii_basic/
    README.md
    cpp/input/model.stl
    cpp/script/
    python/input/model.stl
    python/expected.json
  binary_basic/
    cpp/
    python/
  fixtures/                 # Shared unit-test fixtures
  tools/                    # Shared inspection utilities
  run_case.py               # Python reference-case runner
  stl_to_gmsh.py            # Reusable Python STL geometry implementation
```

- `<case>/cpp` is the native OneFLOW case; its `input/` and `script/`
  belong to that implementation.
- `<case>/python` is the independent Python reference implementation and
  contains `input/model.stl` and `expected.json`.
- Cases may contain only one implementation when a second implementation
  would not add meaningful coverage.
- Generated results are written only to staged runtime workspaces, never to
  the repository source tree.
- Shared fixtures belong in `fixtures/`, not in the case runtime folders.

## Python reference cases

The GUI's **Python reference cases** runner stages each catalogued case and
invokes:

```sh
python grid/stl/run_case.py --case-directory <staged-case-directory>
```

The report is written to the staged case's `results/stl_report.json`.
Expected diagnostic failures pass when the observed diagnostic matches
`expected.json`.

## STL inspection

The reusable reader supports ASCII and binary STL, validates finite
coordinates, and extracts one closed planar boundary loop when supported.
The inspection utility is:

```sh
python grid/stl/tools/inspect.py --input path/to/model.stl --output path/to/report.json
```

STL defines geometry and triangle topology, not CFD boundary conditions.
Boundary names come from an explicit JSON sidecar in the consuming case;
mapping them to CGNS or OneFLOW solver boundary types is a separate task.

## STL to Gmsh

STL-to-Gmsh mesh-generation cases remain a separate workflow under
`grid/stl_to_gmsh`. See that directory's README.

## Tests

```sh
python -m unittest grid.stl.test_stl_to_gmsh
python tools/test_mesh_cases.py validate
```

Tests that generate reports or meshes use temporary directories.

## Scope

Arbitrarily oriented plane projection, multiple disconnected bodies, internal
holes, non-planar surface repair, and a general-purpose production Mesh Core
STL reader are not claimed by this first milestone.
