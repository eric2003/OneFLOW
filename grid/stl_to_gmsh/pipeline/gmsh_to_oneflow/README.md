# Gmsh to OneFLOW pipeline Case

This Case is the Gmsh-to-OneFLOW stage of the STL-to-CFD development workflow.
It consumes the CGNS interchange artifact produced by the preceding Gmsh-to-CGNS
stage, then writes the native OneFLOW grid through the standard Grid task.

The workflow runner supplies sourceGridCaseDir and sourceGridFileName from
the preceding Gmsh-to-CGNS stage. CGNS remains an intermediate artifact rather than the name of this end-to-end stage.
The source CGNS file is therefore not
duplicated into this Case. The generated native grid and grid information
report remain in this Case's results directory.

Run it as part of the workflow:

    python grid/stl_to_gmsh/develop_test/run_cases.py --pipeline-only

The source Case is copied into the user's .oneflow workspace. Its checked-in
parameter files are never rewritten.
