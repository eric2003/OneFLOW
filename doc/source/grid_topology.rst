Grid Topology and Validation
=============================

Purpose and implementation status
---------------------------------

This page describes the topology contracts that should remain explicit as OneFLOW's mesh implementation evolves. It combines the topology-specific material from the repository's engineering note with the broader :doc:`mesh_architecture` overview.

The details below distinguish observed implementation behavior from architectural limitations. In particular, rebuilding derived connectivity is currently an explicit call-path responsibility; it must not be mistaken for automatic cache invalidation.

Topological meaning by dimension
--------------------------------

One dimension
~~~~~~~~~~~~

For a simple chain of two interval cells, there are three nodes, three point-faces, and two cells:

.. code-block:: text

   node:       0           1           2
               o-----------o-----------o
   face:       f0          f1          f2
   cell:           c0          c1

   face       owner cell       neighbour / boundary
   f0         c0               boundary / ghost
   f1         c0               c1
   f2         c1               boundary / ghost

A 1D face is a point, not a line segment. A cell is the interval between nodes. A specialized compact representation is acceptable if it can be converted losslessly to a common face/owner/neighbour interface.

Two dimensions
~~~~~~~~~~~~~~

A 2D mesh requires node, edge/face, and cell connectivity. Shared edges must identify their adjacent cells consistently; boundary edges need a valid boundary classification. Edge-node ordering and owner/neighbour ordering must agree with the normal convention used by numerical fluxes.

Three dimensions
~~~~~~~~~~~~~~~~

For a general polyhedral finite-volume mesh, each face is an ordered loop of nodes, and cells are bounded by faces. A cell's vertex list can support cell-oriented algorithms, but it does not replace oriented face connectivity when surface normals and flux traversal matter. Explicit edge storage is only required when algorithms consume it directly.

Across all dimensions, keep geometric dimension, element type, physical boundary, ghost entity, and partition interface as separate concepts. A single ambiguous integer convention must not represent all of them.

Canonical and derived topology
------------------------------

The intended roles are:

* Node coordinates describe canonical geometry.
* Cell type and cell-to-node connectivity are canonical for cell-centric algorithms and interchange formats.
* Ordered face-to-node connectivity, owner/left cell, neighbour/right cell, and boundary classification are canonical for the current face-based finite-volume model.
* Boundary labels and face-to-region mappings preserve physical meaning; they are not visualization colors.
* Cell-to-face and cell-to-cell adjacency are derived structures when they can be regenerated from canonical connectivity.
* Ghost cells and partition/interface mappings belong to execution and distribution semantics, not to the physical input mesh.

See :doc:`mesh_architecture` for the mesh information levels, snapshot boundary, and format-adapter design.

Boundary normalization contract
-------------------------------

`UnsGrid::NormalizeBc()` converts imported boundary-face neighbour placeholders into local ghost-cell IDs. For a boundary face with index `iFace`, the normalized right-cell ID is:

.. code-block:: text

   ghostCellId = nCells + iFace

The implementation validates real-cell references before changing boundary neighbour IDs. Tests should preserve this failure atomicity: invalid input must not leave only part of the boundary topology normalized.

The formula is an implementation convention, not a portable interchange-format convention. Exported data must document identifier base and scope explicitly.

Derived adjacency and cache refresh
-----------------------------------

The current model has distinct source and derived data:

* `FaceTopo` owns face-to-node and face-to-cell source connectivity plus boundary classification.
* `CellTopo` stores cell-oriented source data and derived `c2f` and `c2c` adjacency arrays.
* `CellTopo::CalcC2f()` and `FaceTopo::CalcC2C()` rebuild derived connectivity from the current source arrays.
* Relevant call paths explicitly rebuild adjacency before reading it, including solver initialization and the METIS partition graph path.

These are explicit rebuild operations, not automatic invalidation. Source containers remain publicly mutable, so editing connectivity and then reading old adjacency without rebuilding it can produce stale results. A revision counter alone would not solve this while callers can mutate arrays through exposed references.

The current targeted rebuild behavior includes:

* `CellMesh::CalcCellSpan()` rebuilds `c2f` and recalculates spans instead of returning early merely because a span array is non-empty.
* `ScalarGrid::CalcC2C(EList&)` reconstructs adjacency from current face connectivity rather than trusting a non-empty caller-provided list.
* `FaceTopo::ReorderLink()` applies one permutation to populated face-indexed arrays, including face types, face-position arrays, and flags, and clears `c2f` because it stores face IDs. Boundary-condition records remain in their original boundary-face order because the permutation preserves the relative order of boundary faces.

These targeted fixes do not establish repository-wide automatic cache coherence. Callers must still rebuild derived data after changing source connectivity. Broader invalidation should follow an audit of actual direct reads and a clear boundary around topology mutation, rather than adding a revision counter that cannot observe direct edits.

Validation invariants
---------------------

Use small hand-constructed meshes with known expected connectivity. At minimum, verify:

#. Counts agree with the arrays they describe; all node and cell references are in range.
#. An internal face connects two distinct real cells.
#. Before ghost normalization, a boundary face references exactly one real cell and has a valid boundary classification.
#. Face-node ordering and owner/neighbour ordering obey the solver's orientation convention.
#. Face-indexed metadata stays aligned after reordering.
#. Interface maps refer to valid faces and zones and use documented ID scopes.
#. A read-normalize-write-read round trip preserves counts, coordinates, face connectivity, internal owner/neighbour pairs, boundary labels, and interface mappings when the format claims to preserve them.
#. Tiny 1D and higher-dimensional fixtures verify adjacency, orientation, positive cell measures where defined, and shared-face traversal suitable for conservation checks.

A useful test ladder is local connectivity tests, boundary-normalization tests, derived-data rebuild tests, reorder consistency tests, round-trip conversion tests, and finally full CFD regressions. Solver output comparisons complement these tests; they do not replace geometry and topology validation.

Known limitation and next audit
-------------------------------

The topology containers still expose mutable source arrays. The implementation therefore relies on explicit rebuild sequencing. Before introducing a general invalidation API, audit direct reads of `c2f`, `c2c`, and related derived values outside the known solver-initialization and partition paths. Any change should be driven by a demonstrated call path and covered by a focused regression test.
