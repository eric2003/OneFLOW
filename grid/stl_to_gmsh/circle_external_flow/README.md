# STL to Gmsh: 2D Circle External Flow (Native C++)

This native OneFLOW Grid Case reads the same 32-segment circular STL as the
Python reference Case. The configured rectangular outer domain encloses the
body; the C++ converter builds an outer curve loop and an oppositely oriented
inner body loop, then creates a planar fluid surface with a hole.

Physical groups are `inflow`, `outflow`, `farfield`, `wall` and `fluid`.
The upper and lower outer boundaries intentionally share the `farfield`
group. These are geometry labels; the Gmsh-to-CGNS stage maps recognized names
to CGNS boundary types, while the solver still owns the actual boundary state.

Run from the repository root after building OneFLOW and installing Gmsh:

    build/bin/Debug/OneFLOW.exe r grid/stl_to_gmsh/circle_external_flow

The generated MSH 4.1 mesh, Gmsh geometry and validation report are written
under this Case's `results/` directory. The mesh is unstructured triangular;
it is not expected to match the Python reference's exact element count.
