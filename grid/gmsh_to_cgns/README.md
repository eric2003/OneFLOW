# Gmsh to CGNS bootstrap Case

This case starts the production OneFLOW grid-conversion workflow from a
Gmsh ASCII MSH 2.x mesh. The initial writer emits a text placeholder with a
`.cgns` suffix. It is deliberately marked `NOT_A_CGNS_FILE`; CGNS tools
must not open or consume it.

The initial reader validates the mesh-format header and counts node and
element records. It does not yet preserve coordinates, element connectivity,
physical groups, or boundary names. Binary MSH and MSH 4.x are rejected.

The Case parameters are in `script/grid.txt`:

```cpp
string sourceGridType = "gmsh";
string targetGridType = "cgns";
```

Copy this Case to an external workspace before running it so generated files
do not pollute the source tree:

```text
OneFLOW.exe r <workspace>/gmsh_to_cgns
```

Next steps: introduce a format-neutral mesh model, map Gmsh physical groups
to named boundaries, then implement a real CGNS writer and CGNS-reader
round-trip tests. No Gmsh API is used in this first stage.
