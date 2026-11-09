# SLAU2 CFD pipeline Case

This dedicated Solve Case consumes the native OneFLOW grid produced by the
preceding CGNS-to-OneFLOW stage. The workflow runner injects sourceGridCaseDir
and rewrites gridFileName in the staged copy of script/cfd.txt; the source
parameter files are not modified.

The configuration uses the supplied 50-step SLAU2 parameter set, including
global_dt = 0.001 and LU-SGS time integration. It deliberately sets idump = 0
instead of the flat-plate template's idump = 1: this square has inflow,
outflow, symmetry and farfield boundaries, but no wall on which a laminar
flat-plate initial profile would be based. Residual, restart and visualization
outputs remain enabled. The pipeline requires non-empty:

- results/res.dat
- restart/flow.dat
- visual/visual.dat

The run is an end-to-end runtime smoke test. It does not by itself establish
numerical accuracy, convergence, conservation or grid independence.
