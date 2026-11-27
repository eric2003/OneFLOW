OECD/NEA–Vattenfall T-junction Case
====================================

This Case reproduces the three-dimensional fluid-domain geometry and published
operating point of the OECD/NEA–Vattenfall T-junction benchmark.

Benchmark definition
--------------------

The test section uses a 140 mm main-pipe inner diameter and a 100 mm branch
inner diameter. The main straight-pipe extensions beyond the junction block are each 1070 mm.
Figure A2 gives 1220 mm from the junction midpoint to each main pipe end
(1070 mm straight-pipe length plus the 150 mm centreline-to-extension
reference shown in that drawing). Table A1 separately records the Plexiglas
block as 220 x 285 x 325 mm. The branch straight-pipe
extension above the block is 470 mm; the block top is 223 mm above the main
centreline, giving 693 mm from the centreline to the branch inlet. The baseline
fluid domain is the Boolean union of the cylindrical bores. The intersection
edge was described as slightly rounded, with a radius below 0.2 mm that was
not measured. The baseline therefore uses a sharp intersection and does not
invent a fillet.

The benchmark operating point is 9.0 L/s at 19 °C in the main inlet and 6.0 L/s
at 36 °C in the branch inlet. The temperature difference is 17 K. The
benchmark report specifies de-ionized tap water and points to the 1967 ASME
Steam Tables for properties, with thermal-expansion data from Incropera and
DeWitt.

Geometry, parameters and equations
----------------------------------

The machine-readable source of truth is
`grid/stl_to_gmsh/3d/pipe_tee/oecd_nea_vattenfall/automatic_tetrahedra/benchmark.json`.
It records dimensions, operating conditions, derived bulk velocities and
Reynolds-number estimates (using linearly interpolated Table A2 properties),
water/Plexiglas properties, measurement locations, uncertainties, formulas and
the future CFD-domain requirements.

The principal relations are

.. math::

   A = \frac{\pi D^2}{4}, \qquad
   U_b = \frac{Q}{A}, \qquad
   Re_D = \frac{\rho U_b D}{\mu}

and

.. math::

   T^* = \frac{T-T_{cold}}{T_{hot}-T_{cold}}.

The main inlet profile is fully developed and measured at x=-3D. The branch
profile is developing and measured 3.1 branch diameters upstream of the
junction. The report notes that the supplied inlet profiles came from an
earlier test with a different main/branch flow-rate ratio; this provenance is
retained in the benchmark data contract.

Reference measurements
----------------------

The report documents upstream LDV profiles at -3D_main and -3.1D_branch,
downstream LDV stations at 2.6D and 6.6D, downstream PIV velocity measurements
at 1.6D, 2.6D, 3.6D and 4.6D, and wall thermocouple measurements from 2D to 20D
downstream. The PIV measurements are isothermal; thermocouple data characterize
the thermal mixing case. The report also defines a curve-comparison metric
based on the mean absolute difference between computed and measured values.

Experimental data ownership
---------------------------

The original experimental data are owned by Vattenfall Research and
Development AB. The report states that tabular files were distributed to
benchmark participants. OneFLOW records their existence and provenance but
does not redistribute those files until source files and permission are
confirmed.

Current OneFLOW scope
---------------------

The Python reference and C++17 geometry driver use Gmsh/OpenCASCADE to create
one fluid-domain BooleanUnion from the perpendicular main and branch bores.
The source surface mesh keeps `wall`, `main_inflow`, `branch_inflow` and
`outflow` distinct, and the geometry exporter preserves checksum-bound STL
patches as portable surface evidence. The final MSH4 tetrahedral mesh is
generated directly from the same OpenCASCADE geometry script, not rebuilt from
the STL patches. This preserves both inlet Physical Groups and avoids a second
discrete-surface topology reconstruction. The shared MSH4 inspector checks
every tetrahedron vertex, edge midpoint and centroid against the analytic
union of the two cylindrical bores, in addition to checking boundary
connectivity and Physical Groups. A checksum-bound `*.geometry-validation.json`
report is written beside the final MSH4. Mesh Project creation, replacement and
import preserve the T-junction Case identity and inlet subregions. The Case is
registered in the existing Case Catalog and Mesh-to-CFD Validation workflow.

The current CFD stage is a bounded two-step smoke test only. The generated mesh is the local test section: its inlet caps are at the pipe ends, not at the benchmark LDV measurement planes, and its outlet does not reach the 20D thermocouple station. The future CFD-domain crop/extension requirements are explicitly recorded in benchmark.json. It does not establish
mesh convergence, physical accuracy, thermal mixing fidelity, turbulence
statistics or agreement with the benchmark. Future physical validation must
record fluid properties, inlet profiles, turbulence and thermal models, wall
treatment, mesh refinement, time step, averaging duration and data uncertainty.

References
----------

* B. L. Smith, J. H. Mahaffy and K. Angele, *Report of the OECD/NEA-Vattenfall
  T-Junction Benchmark exercise*, NEA/CSNI/R(2011)5 (2011).
  https://www.oecd-nea.org/upload/docs/application/pdf/2021-03/csni-r2011-5.pdf
* B. L. Smith, J. H. Mahaffy and K. Angele, “A CFD benchmarking exercise based
  on flow mixing in a T-junction,” *Nuclear Engineering and Design*, 264,
  80–88 (2013). DOI: 10.1016/j.nucengdes.2013.02.030.
