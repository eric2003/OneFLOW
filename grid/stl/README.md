# STL format I/O and compatibility cases

This directory covers STL file handling and focused compatibility cases. It is
not the STL-to-Gmsh mesh-generation workflow; that workflow and its native C++
and Python Cases live in `grid/stl_to_gmsh/`.

## Case layout

Each compatibility Case is organized by the input property or behavior it
checks. Native C++ and Python reference implementations are separated when
both provide useful, independent evidence.

```text
grid/stl/
  ascii_basic/
    cpp/input/model.stl
    cpp/script/
    python/input/model.stl
    python/expected.json
  binary_basic/
    cpp/
    python/
  fixtures/                 # Shared inputs for focused Python unit tests
  inspect.py                # Standalone STL inspection command
  run_case.py               # Python reference-case runner
  stl_to_gmsh.py            # Shared planar geometry/Gmsh helper used by conversion Cases
  test_stl_to_gmsh.py       # Unit tests for the Python geometry helper
```

- `<case>/cpp` is a native OneFLOW Grid Case. Its `input/` and `script/`
  belong to that Case.
- `<case>/python` is an independent Python reference Case. Its input and
  expected diagnostics are local to that Case.
- A second implementation is added only when it tests a distinct path or
  provides a useful independent cross-check.
- Generated results belong in an external staged workspace, never in the
  repository source tree.

## Shared test fixtures

`fixtures/` is deliberately limited to small, stable inputs shared by the
focused unit tests in `test_stl_to_gmsh.py`. For example, `square.stl` and
`square.boundary.json` form a minimal planar-geometry test pair.

These are not runtime Case inputs. Runtime Cases keep their own
`<case>/<implementation>/input/` files so each Case remains independently
runnable and understandable. Do not add a fixture here merely to avoid
copying a file used by only one Case.

## STL inspection

`inspect.py` is a command-line diagnostic utility. It reads an explicitly
provided STL and writes a JSON report to an explicitly provided destination:

```sh
python grid/stl/inspect.py --input path/to/model.stl --output path/to/report.json
```

The report includes file format, triangle count, coordinate bounds, XY
planarity, and whether the supported single closed planar boundary can be
reconstructed. The caller chooses the output path; the tool must not create
runtime output under the source tree by default.

## Python reference cases

The GUI's **Python reference cases** runner stages each catalogued Case and
invokes:

```sh
python grid/stl/run_case.py --case-directory <staged-case-directory>
```

The report is written to the staged Case's `results/stl_report.json`.
Expected diagnostic failures pass when the observed diagnostic matches
`expected.json`.

## STL-to-Gmsh conversion

The reusable Python geometry helper in this directory is consumed by the
separate conversion Cases under `grid/stl_to_gmsh/`. It is an implementation
dependency, not a second set of conversion Cases. See
`grid/stl_to_gmsh/README.md` for the official progression and execution
instructions.

## Tests

```sh
python -m unittest grid.stl.test_stl_to_gmsh
python tools/test_mesh_cases.py validate
```

Tests that generate reports or meshes use temporary directories.

## Scope

The compatibility Cases establish specific STL file and geometry behaviors.
They do not claim arbitrary plane projection, multiple disconnected bodies,
internal-hole repair, non-planar surface repair, or a general-purpose
production Mesh Core topology service.
