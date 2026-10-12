# Grid topology model and validation plan

This note records a working model for OneFLOW's grid architecture. It is a guide for incremental refactoring, not a proposal to replace all current data structures at once.

## 1. The minimum topological meaning

### One-dimensional chain

For a simple chain of two cells, there are three nodes, three faces (points in 1D), and two cells:

```text
node:       0           1           2
            o-----------o-----------o
face:       f0          f1          f2
cell:           c0          c1

face       left/owner cell       right/neighbour cell
f0         c0                     boundary / ghost
f1         c0                     c1
f2         c1                     boundary / ghost
```

The geometry of a 1D face is a point, not a line segment. The cell is the interval between two nodes. For a simple ordered chain, cell-to-face relationships can be inferred from ordering, but a general solver interface can still expose the same owner/neighbour contract as higher dimensions.

### Two-dimensional example

Two cells sharing one edge:

```text
n2 ------------- n3 ------------- n4
|                |                |
|    cell 0      |     cell 1      |
|                |                |
n0 ------------- n1 ------------- n5
                 |
              f_shared
```

The sketch is schematic: each cell has an ordered node loop, each edge has an ordered node list, and the shared edge has two adjacent cells. Boundary edges have one real adjacent cell and a boundary label/patch.

### Three-dimensional example

For a polyhedral finite-volume mesh, a face is an ordered loop of nodes; a cell is bounded by faces. The face-to-cell relation is normally represented by an owner and an optional neighbour. A cell's full vertex list may be useful for other algorithms, but it is not a substitute for oriented face connectivity when fluxes and surface normals matter.

## 2. What should be canonical versus derived?

| Information | Proposed role | Reason |
| --- | --- | --- |
| Node coordinates | Canonical geometry | Needed for metrics and geometric queries |
| Cell type and cell-to-node connectivity | Canonical for cell-centric formats and element algorithms | Required by many import/export formats |
| Face-to-node ordered connectivity | Canonical for the current finite-volume face model | Preserves face orientation and supports normals/fluxes |
| Face owner/left cell and neighbour/right cell | Canonical for face-based finite-volume traversal | Makes internal versus boundary face semantics explicit |
| Boundary labels and face-to-patch mapping | Canonical boundary metadata | Connects geometric boundary entities to physical boundary conditions |
| Cell-to-face / cell-to-cell adjacency | Derived where practical | Avoids duplicated mutable truth; cache only when measured performance justifies it |
| Ghost cells and partition/interface maps | Execution/distribution layer | These are not additional physical cells in the input mesh |
| Face centers, normals, cell volumes, search trees | Derived geometry/cache | Rebuild or invalidate when geometry/topology changes |

This is intentionally a hybrid recommendation. OpenFOAM uses ordered face points with owner/neighbour cell labels and boundary patches. MFEM provides element-to-vertex connectivity and tables for element/face/edge relations. VTK uses cell types plus connectivity and offsets for general mixed-cell representation. These are different designs serving different algorithms; OneFLOW should borrow their clear invariants, not copy one implementation wholesale.

## 3. Dimension does not dictate storage layout

- **1D:** node, interval-cell, and point-face semantics are enough for a simple line mesh. A specialized compact representation is valid if it can convert losslessly to the common face/owner/neighbour interface.
- **2D:** node, edge/face, and cell connectivity are needed; edge orientation and boundary classification matter.
- **3D:** node, face, and cell connectivity are needed for a general polyhedral finite-volume mesh. Explicit edge storage is optional unless algorithms consume it directly.
- **All dimensions:** distinguish geometric dimension, element type, physical boundary, ghost entity, and partition interface. Do not encode these concepts through one ambiguous integer convention.

## 4. Invariants and test ladder

1. **Local connectivity:** counts agree; all node/cell IDs are in range; internal faces connect two distinct real cells; boundary faces connect exactly one real cell before ghost normalization.
2. **Orientation:** a face's node ordering and owner/neighbour ordering agree with the solver's normal convention.
3. **Boundary metadata:** each boundary face has a valid boundary classification; interface maps point to valid boundary faces and zones.
4. **Round-trip conversion:** read -> normalize -> write -> read preserves counts, coordinates, face-node connectivity, internal owner/neighbour pairs, boundary labels, and interface mappings.
5. **Partition semantics:** each local interface face maps to the intended remote zone/face; ghost IDs are not confused with global physical cell IDs.
6. **Physics checks:** for tiny 1D and degenerate higher-dimensional meshes, verify expected cell adjacency, face orientation, positive cell measures, and conservation-compatible shared-face traversal.

The unit tests should use tiny hand-constructed meshes with known answers. File-based CFD regressions remain valuable, but a solver output comparison alone is not a complete geometry/topology test.

## 5. Immediate scope

The current `UnsGrid::NormalizeBc()` contract is tested directly: validate all real-cell references before mutating boundary neighbour IDs, then assign boundary ghost IDs as `nCells + boundaryFaceIndex`. The tests intentionally check that a failed validation does not partially normalize boundary faces.

The next architecture step should audit the conventions at each importer/converter boundary and build a read-back topology validator before changing ownership or replacing storage containers.
