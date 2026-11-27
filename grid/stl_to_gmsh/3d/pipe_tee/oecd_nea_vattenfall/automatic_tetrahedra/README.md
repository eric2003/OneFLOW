# OECD/NEA–Vattenfall T-junction: automatic tetrahedral mesh

This Case adds a literature-backed three-dimensional T-junction to the existing
STL-to-Gmsh and Mesh-to-CFD workflow. It follows the same staged Case pattern as
the 90-degree pipe bend; it does not introduce a second meshing or validation
framework.

Case ID: `stl_to_gmsh_pipe_tee_oecd_nea_vattenfall`

## Benchmark identity and geometry

The geometry and operating point are based on the OECD/NEA–Vattenfall T-junction
benchmark (NEA/CSNI/R(2011)5, 2011):

- Main pipe inner diameter: 140 mm.
- Branch pipe inner diameter: 100 mm.
- Main straight-pipe extensions beyond the junction block: 1070 mm each.
- Figure A2 gives 1220 mm from the junction midpoint to each main pipe end
  (1070 mm straight-pipe length plus the 150 mm centreline-to-extension
  reference shown in that drawing). Table A1 separately records the Plexiglas
  block dimensions as 220 x 285 x 325 mm.
- Branch straight-pipe extension above the junction block: 470 mm. Figure A2
  places the block top 223 mm above the main centreline, giving 693 mm from the
  centreline to the branch inlet.
- Main inlet: 9.0 L/s at 19 °C.
- Branch inlet: 6.0 L/s at 36 °C.
- Temperature difference: 17 K.
- The intersection edge is reported as very slightly rounded, with an
  unmeasured radius below 0.2 mm. The baseline geometry deliberately uses a
  sharp intersection rather than inventing a fillet radius.

The fluid domain is the Boolean union of the three cylindrical bores. The main
pipe axis is +X, the branch inlet axis is +Z, and the origin is at the
centreline intersection. The source surface MSH2 preserves dimension-2 `wall`, `main_inflow`,
`branch_inflow` and `outflow` as distinct groups. The shared volume mesher
currently combines the two inlet patches into the compatibility `inflow`
region, so the final MSH4 contains `wall`, `inflow`, `outflow` and the
`fluid_domain` volume group. The separate checksum-bound inlet STL subregions
are retained for future benchmark-specific boundary conditions.

The initial geometry producer uses Gmsh/OpenCASCADE to form the union, mesh its
boundary, and export a closed STL plus checksum-bound region patch files. The
existing generic region-aware STL-to-Gmsh mesher then generates the final
TETRA_4 volume mesh. The Gmsh surface-mesh step is a geometry discretization
step; it is not a substitute for the shared final volume-mesh stage.

## Current delivery scope

- Reuse the shared STL region manifest, generic volume mesher, MSH4 inspector,
  Gmsh-to-CGNS converter, CGNS-to-OneFLOW converter and CFD smoke consumer.
- Preserve wall, both inlet patches (combined as `inflow`) and the main outlet
  as explicit semantic regions.
- Run only a bounded CFD smoke test through the existing Mesh-to-CFD
  Validation Workbench. The current smoke configuration is not a benchmark
  calculation.
- Keep generated STL, intermediate surface MSH, final volume MSH, reports and
  validation artifacts in the external Case workspace, never in the source tree.

## Reference data and future physical validation

The machine-readable reference contract is
[`benchmark.json`](benchmark.json). It records the dimensions, flow rates,
temperatures, derived bulk velocities, estimated Reynolds numbers derived from linearly interpolated Table A2 properties, formulas,
measurement stations, reference quantities, source documents and future
comparison plan.

Important benchmark details:

- The main inlet profile is fully developed and measured at x = -3D.
- The branch inlet profile is developing and measured 3.1 branch diameters upstream (z = +0.31 m in the OneFLOW coordinate convention).
- The published benchmark inlet profiles came from an earlier flow-ratio-2
  test; the branch flow remained 6 L/s while the main flow was 9 L/s for the
  benchmark test. This detail must not be lost when future inlet profiles are
  reconstructed.
- Downstream LDV stations are at 2.6D and 6.6D. Downstream PIV velocity
  measurements are isothermal at 1.6D, 2.6D, 3.6D and 4.6D. Thermocouple data
  cover 2D through 20D downstream.
- The report describes a comparison metric
  `M = (1/N) sum_i |C_i - D_i|`.

The benchmark report states that the working fluid is de-ionized tap water and
references the 1967 ASME Steam Tables for properties, with thermal expansion
data from Incropera and DeWitt. The actual selected property values and
interpolation states must be recorded in a future physical-validation run.

The experimental data are owned by Vattenfall Research and Development AB.
The report documents the measurement quantities and locations and notes that
tabular data files were distributed to benchmark participants. Original
tabular data are not copied into this repository until their source and
redistribution permissions are confirmed. The report and journal reference are
recorded in `benchmark.json`.

## Verification status

The lifecycle contract deliberately separates implementation from evidence:

- Geometry generation and named surface regions: implementation target.
- Final MSH4 tetrahedral mesh and Physical Groups: implementation target.
- CGNS conversion, OneFLOW ingestion and bounded CFD smoke: shared workflow
  stages, to be marked verified only after actual run evidence.
- Mesh-quality acceptance and physical validation: not completed.
- Experimental comparison, thermal statistics and spectra: reserved for
  future work.

A successful short smoke run does not demonstrate convergence, correct
turbulence statistics, thermal mixing accuracy, or agreement with the
Vattenfall experiment.

## References

1. B. L. Smith, J. H. Mahaffy and K. Angele, *Report of the OECD/NEA-Vattenfall
   T-Junction Benchmark exercise*, NEA/CSNI/R(2011)5, 2011.
   https://www.oecd-nea.org/upload/docs/application/pdf/2021-03/csni-r2011-5.pdf
2. B. L. Smith, J. H. Mahaffy and K. Angele, “A CFD benchmarking exercise based
   on flow mixing in a T-junction,” *Nuclear Engineering and Design*, 264,
   80–88 (2013). DOI: 10.1016/j.nucengdes.2013.02.030.


### Preserved inlet subregions

The source surface mesh has separate `main_inflow` and `branch_inflow`
Physical Surfaces. Before exporting the STL, the shared Python exporter
validates the three port planes, circular areas, area centroids, and wall-node
distances against the two perpendicular cylinder axes. The measured results
are stored in the `geometry_validation` section of the geometry manifest. The
native C++ driver passes its geometry parameters through the same validator,
so both generation paths enforce the same geometric contract. The generator
also writes separate checksum-bound `main_inflow.stl` and
`branch_inflow.stl` patches. The shared volume-mesher
contract currently accepts the three-region `wall`/`inflow`/`outflow`
partition, so it combines the two inlet patches into the compatibility
`inflow` group for the current smoke-test MSH4. Their separate source patches
remain available for future benchmark-specific inlet profiles and temperatures.
