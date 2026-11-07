# Two Adjacent Quadrilaterals

This is the first self-contained mesh in the topology-oriented grid
conversion Case collection. It exercises a two-cell interface rather than
only proving that a single element can be read.

## Geometry and topology

The MSH 4.1 ASCII input describes a unit-height, two-unit-wide rectangle split
at x = 1.

| Quantity | Expected |
|---|---:|
| Topological dimension | 2 |
| Nodes | 6 |
| Top-dimensional cells | 2 QUAD_4 |
| Geometric curve entities | 7 |
| Boundary LINE_2 elements | 6 |
| Internal interface curve / line element | 1 |
| External physical boundary groups | 4 |
| Total Gmsh elements | 9 |

The two quadrilaterals are (1, 2, 5, 4) and (2, 3, 6, 5). They share the edge
between nodes 2 and 5. The shared edge is represented as a geometric curve
and a Gmsh line element, but it is deliberately not assigned to an external
boundary Physical Group.

| Physical Group | Edges | Intended meaning |
|---|---|---|
| inflow | Left | Inlet boundary |
| outflow | Right | Outlet boundary |
| symmetry | Bottom, split into two edges | Symmetry boundary |
| farfield | Top, split into two edges | Farfield boundary |
| fluid_domain | Both quadrilateral surfaces | 2-D fluid region |

The converter must use Physical Group names, not infer boundary meaning from
edge coordinates. The shared interface must not be exported as an external
boundary condition.

## Conversion route

The checked-in script/grid.txt selects Gmsh MSH 4.1 ASCII to CGNS.
The Case writes results/two_quads.cgns and the mesh information report
configured in the same parameter file. These are runtime artifacts and are
not checked into the repository.

## Running

1. Build OneFLOW with CGNS support.
2. Start the desktop GUI from the repository root: python -m tools.oneflow_gui.
3. Select Gmsh Two Adjacent Quadrilaterals in the Case Workbench and run it.
   The GUI executes a protected workspace copy; the checked-in input remains
   unchanged.
4. Inspect the generated CGNS file and conversion/read reports in the
   execution workspace.

The catalog preflight can be checked separately:

    python tools/mesh_cases.py validate two_quads

This command validates Case registration and parameter-file references only;
it does not execute the conversion. The Gmsh-to-CGNS example remains a manual
GUI workflow and is not part of CTest, Python unittest, or GitHub Actions.

## What this example proves

- A MSH 4.1 file can contain more than one top-dimensional quadrilateral.
- Two cells can share an internal edge while the exterior boundary remains
  explicitly classified.
- One physical group may cover multiple geometric curve entities.
- Boundary semantics remain separate from cell adjacency.

## What it does not prove

This planar example does not establish support for mixed element families,
high-order elements, periodicity, partition/ghost metadata, or 3-D volume
cells. Those require separate runnable Cases with their own inputs and
expected topology.
