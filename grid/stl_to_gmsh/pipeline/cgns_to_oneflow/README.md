# CGNS to OneFLOW pipeline Case

This Case is the CGNS consumer stage of the STL-to-CFD development workflow.
It uses the standard OneFLOW Grid task and the existing cgns-to-oneflow
conversion route.

The workflow runner supplies sourceGridCaseDir and sourceGridFileName from
the preceding Gmsh-to-CGNS Case. The source CGNS file is therefore not
duplicated into this Case. The generated native grid and grid information
report remain in this Case's results directory.

Run it as part of the workflow:

    python grid/stl_to_gmsh/develop_test/run_cases.py --pipeline-only

The source Case is copied into the user's .oneflow workspace. Its checked-in
parameter files are never rewritten.
