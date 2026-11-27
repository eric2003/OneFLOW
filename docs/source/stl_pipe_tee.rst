OECD/NEA–Vattenfall T-junction Case
====================================

This Case reproduces the three-dimensional fluid-domain geometry and published
operating point of the OECD/NEA–Vattenfall T-junction benchmark.

Benchmark definition
--------------------

The test section uses a 140 mm main-pipe inner diameter and a 100 mm branch
inner diameter. The main straight-pipe extensions beyond the junction block are each 1070 mm.
Figure 2 places each block side 150 mm from the main-pipe centreline, giving
1220 mm from the centreline to each main pipe end. The branch straight-pipe
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
It records dimensions, operating conditions, derived bulk velocities, reported
Reynolds numbers, measurement locations, formulas and the future validation
plan.

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
the fluid-domain Boolean union, then export the same closed STL and named
surface patches. The shared region-aware STL mesher creates the final MSH4
tetrahedral mesh. The Case is registered in the existing Case Catalog and
Mesh-to-CFD Validation workflow.

The current CFD stage is a bounded smoke test only. It does not establish
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
