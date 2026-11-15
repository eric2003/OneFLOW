# Gmsh 2D to OneFLOW Case

This case demonstrates the supported two-dimensional conversion route:

`Gmsh MSH 4.1 ASCII -> CGNS -> OneFLOW unstructured grid`.

The route deliberately uses CGNS as an independently readable interchange
artifact. It is not a separate direct Gmsh-to-runtime mesh implementation.
The same MSH reader and CGNS writer are therefore exercised by both
`gmsh -> cgns` and `gmsh -> oneflow`.

## Mesh contract

The input `grid/single_quad.msh` is a minimal one-cell quadrilateral:

| Item | Expected |
|---|---:|
| Topological dimension | 2 |
| Nodes | 4 |
| Volume cells | 1 QUAD_4 |
| Boundary elements | 4 BAR_2 |
| Physical boundary groups | 4 |

The four curve Physical Groups are explicit semantic labels:

| Gmsh Physical Group | CGNS boundary type | OneFLOW boundary type |
|---|---|---|
| `inflow` | `BCInflow` | `INFLOW` |
| `outflow` | `BCOutflow` | `OUTFLOW` |
| `symmetry` | `BCSymmetryPlane` | `SYMMETRY` |
| `farfield` | `BCFarfield` | `FARFIELD` |

Names are read from Gmsh Physical Names; the converter does not infer boundary
meaning from edge coordinates. The mesh is intentionally small so that all
four faces are boundaries and the conversion can be inspected directly.

## Preflight

From the repository root:

```sh
python grid/mesh_workflow/tools/validate_case.py grid/gmsh_to_oneflow
```

The validator checks the case control-file includes, required grid parameters,
the Gmsh input path, supported source/target format pair, and output
directories. It does not replace OneFLOW runtime validation.

## Run without modifying the repository case

Run a workspace copy so generated files never enter the source checkout.

PowerShell:

```powershell
$case = Join-Path $env:USERPROFILE ".oneflow/cases/gmsh_to_oneflow"
New-Item -ItemType Directory -Force (Split-Path $case) | Out-Null
Copy-Item -Recurse -Force ".\grid\gmsh_to_oneflow" $case

# Run from the OneFLOW executable directory so shared system resources resolve.
Push-Location ".\build\bin\Debug"
.\OneFLOW.exe r $case
Pop-Location
```

Linux (adjust the executable path for the selected build configuration):

```sh
case_dir="$HOME/.oneflow/cases/gmsh_to_oneflow"
mkdir -p "$(dirname "$case_dir")"
cp -R grid/gmsh_to_oneflow "$case_dir"
(
  cd build/bin
  ./OneFLOW r "$case_dir"
)
```

The case uses the standard `OneFLOW.exe r <case-directory>` entry point.
The executable working directory should remain its build/runtime directory;
case-relative input and generated outputs are resolved from the case directory.

## Expected outputs

After a successful conversion:

- `results/single_quad.gmsh-to-oneflow.cgns`: CGNS interchange file.
- `results/single_quad.oneflow.grid`: native OneFLOW unstructured grid.
- `results/single_quad.gmsh-to-oneflow.gmsh-read-report.json`: machine-readable
  Gmsh read report.
- `results/single_quad.gmsh-to-oneflow.gmsh-report.md`: human-readable
  inspection report.
- `grid/bcRegionMap.txt`: generated OneFLOW boundary-name mapping.

The CGNS and native OneFLOW files are separate outputs. The native grid writer
must not overwrite the CGNS interchange file.

## Runtime smoke test

The CTest
`GmshSingleCellRuntimeIntegration.ConvertsReadsComputesAndWritesFlowField`
uses this exact one-cell Gmsh input, performs the same Gmsh-to-CGNS-to-OneFLOW
conversion, then runs one production solver iteration with the standard
solver-case configuration. It checks that residual, restart-flow and
visualization outputs are produced.

That test is a runtime workflow smoke test, not a numerical-accuracy claim.
The next numerical-verification stage should use a multi-cell mesh and a
manufactured or analytical solution; this one-cell fixture first protects
mesh ingestion, boundary identity and runtime output lifecycle.
