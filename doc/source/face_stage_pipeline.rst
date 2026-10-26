Face Stage Pipeline
===================

Purpose and status
------------------

This page records the **staged face** construction path used when building
unstructured calculation grids from scanned elements and CGNS boundary data.
It complements :doc:`grid_topology` (runtime topology contracts) and
:doc:`mesh_architecture` (information levels). Implementation and tests take
precedence if this page lags the code.

The page is maintained next to the code changes it describes so the design
handbook grows with the product rather than as a separate late write-up.

Goals
-----

#. One controlled write path for new faces and boundary-condition slots.
#. Parallel face arrays stay length-aligned (topology, flags, BC, child lists).
#. Boundary BC export and face reordering keep a checkable contract.
#. Face *identity* (deduplication) stays separate from face *orientation*
   (node ring and left/right ownership).

Pipeline overview
-----------------

Orchestration is in ``GridElem``; mutation rules live in ``FaceSolver`` and
``FaceTopo``.

::

    Scan volume / polygon / polyhedron elements
        -> FaceSolver::StageNewFace / Claim | Attach | Bind
    Scan BC (CGNS BOCO)
        -> AssignBcToMatchingFaces / SetStagedFaceBc
    GenerateCalcElement
        -> AppendBoundaryBcToRecord
    GenerateCalcGrid
        -> ReleaseFaceTopo
        -> BcTypeMap::MapFromCgnsInPlace
        -> FaceTopo::ReorderLink
        -> ValidateBoundaryBcLayout

Staged face invariant
---------------------

``StageNewFace`` grows the following together:

* ``FaceTopo`` faces, types, flags, left/right cells, and face-local positions
* staged BC slots (private; read through ``GetStagedBcType`` /
  ``GetStagedBcNameId``)
* ``childFid`` slots (private; read through ``GetChildFaceIds``)

``StagedBcMatchesFaceCount`` is the named alignment check.

Post-stage writes:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - API
     - Role
   * - ``SetStagedFaceBc``
     - Sole BC write path after stage
   * - ``ClaimStagedFaceAsLeft`` / ``AttachStagedFaceAsRight``
     - Polygon ownership
   * - ``BindCellToExistingFace``
     - Traditional element faces

**Tests:** ``FaceSolverStageTest``, ``FaceSolverLifecycleTest``.

Boundary BC export and reorder
------------------------------

#. ``AppendBoundaryBcToRecord`` walks faces in current order and, for each face
   with ``right == INVALID_INDEX``, appends staged BC into ``BcRecord``.
#. ``ReorderLink`` applies ``BuildBoundaryFirstPermutation(rightCells)`` so
   boundary faces form a prefix while preserving relative order within the
   boundary group and within the interior group.
#. ``ValidateBoundaryBcLayout`` requires:

   * ``BcRecord`` size equals the boundary face count
   * faces ``[0, nB)`` are boundary and ``[nB, n)`` are interior

Because the permutation preserves boundary relative order, ``BcRecord`` need
not be permuted separately after a correct ``ReorderLink``.

**Tests:** ``FaceTopoReorderLinkTest`` (permutation, layout validation).

Face identity versus orientation
--------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Mechanism
     - Meaning
   * - ``HXKey`` / ``HXLookup`` / ``FaceSearchBasic``
     - Undirected identity: node ids are sorted for deduplication
   * - ``GetFaceNodeRing`` / left-right cells
     - Orientation: first-seen node ring and cell ownership

The same node set in any order must not create a second internal face record.
Flux orientation uses ring order and left/right (or signed references), not a
second identity entry.

**Tests:** ``HXLookupTest``, ``FaceSearchBasicTest``.

BC classification and CGNS mapping
----------------------------------

* After export, CGNS type ids are converted with
  ``BcTypeMap::MapFromCgnsInPlace`` (**test:** ``BcTypeMapTest``).
* ``BC::IsInterfaceBc`` treats ``BC::INTERFACE`` and any negative code as
  interface (historical convention; **test:** ``BcClassificationTest``).

Region names use injectable ``RegionNameRegistry``; ``RegionNameMap`` and
``VolumeNameMap`` provide separate process-wide defaults.
``RegisterRegionNames`` encodes ordered name-id assignment for BOCO-like loops
(**test:** ``RegionNameRegistryTest``).

Interface split helpers
-----------------------

* ``GetInterfaceChildFaceIds(l2g, splitTable, zone, local)`` is the pure
  overload (**test:** ``InterfaceChildFaceIdsTest``).
* ``ValidateGlobalFaceMapping`` / ``ValidateGlobalFaceMappings`` check
  ``gI2Zid`` / ``g2l`` alignment (**test:** ``GlobalFaceMappingTest``).
* ``RegisterLocalInterfaceFace`` records ``l2g`` / ``gI2Zid`` / ``g2l`` when a
  local interface face joins a global id (**test:**
  ``RegisterLocalInterfaceFaceTest``).
* ``ExpandBcRecordForSplitInterfaces`` repeats BC slots when an interface face
  splits into children (**test:** ``ExpandBcRecordSplitTest``).
* ``BuildInverseInterfaceMapForZone`` fills ``gI2ZidNew`` / ``g2lNew`` from a
  zone's ``l2gNew`` row (at most two zone slots) (**test:**
  ``BuildInverseInterfaceMapTest``).
* Quad diagonal candidates:
  ``BuildQuadDiagonalLocalTris`` / ``MapLocalTrisToGlobalNodes``
  (**test:** ``QuadDiagonalSplitTest``).

Open design pressures
---------------------

These remain harder to unit-test and mark follow-up work:

#. ``IFaceLink`` still exposes large public tables (``l2g``, ``gI2Zid``, ...);
   registration is covered by ``RegisterLocalInterfaceFace``, but mutation
   still happens through public fields.
#. Line mid-point split still needs ``PointLocator`` and ``nodeMesh`` (geometry
   services are not yet injectable).
#. End-to-end CGNS zone scan still requires heavyweight zone setup.

When those paths are split into pure helpers or ports, update this page in the
same change and add tests.

PR checklist (executable norm)
------------------------------

Use this list when a change touches staged faces, BC export, or boundary-first
reorder. The goal is one pull request that keeps **code**, **tests**, and
**this page** aligned.

Always
^^^^^^

* Prefer pure helpers or named APIs over new ad-hoc writes to parallel arrays.
* Do not introduce non-ASCII characters in source or in this page.
* If behavior is intentional, name the unit test that locks it; if no test
  exists yet, add one in the same change when the logic can be isolated.
* If a path stays hard to test, add or extend an item under `Open design
  pressures` instead of only describing it in the PR text.

Change type matrix
^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - If you change...
     - Update tests (minimum)
     - Update docs (minimum)
   * - ``StageNewFace`` / parallel growth /
       ``StagedBcMatchesFaceCount``
     - ``FaceSolverStageTest``,
       ``FaceSolverLifecycleTest``
     - `Staged face invariant` on this page
   * - ``SetStagedFaceBc`` /
       ``AssignBcToMatchingFaces`` /
       BC slot privacy or accessors
     - ``FaceSolverStageTest`` (BC write/read)
     - write-path table on this page
   * - ``ClaimStagedFaceAsLeft`` /
       ``AttachStagedFaceAsRight`` /
       ``BindCellToExistingFace``
     - ``FaceSolverStageTest`` (ownership)
     - post-stage writes table on this page
   * - ``AppendBoundaryBcToRecord``
     - stage/export tests that cover boundary
       faces and alignment failures
     - `Boundary BC export and reorder`
   * - ``BuildBoundaryFirstPermutation`` /
       ``ReorderLink`` /
       ``ValidateBoundaryBcLayout``
     - ``FaceTopoReorderLinkTest``
     - `Boundary BC export and reorder`
   * - ``HXKey`` / ``HXLookup`` /
       ``FaceSearchBasic`` identity rules
     - ``HXLookupTest``,
       ``FaceSearchBasicTest``
     - `Face identity versus orientation`
   * - ``BcTypeMap`` or ``BC::Is*`` predicates
     - ``BcTypeMapTest``,
       ``BcClassificationTest``
     - `BC classification and CGNS mapping`
   * - ``RegionNameRegistry`` /
       ``RegisterRegionNames``
     - ``RegionNameRegistryTest``
     - region-name paragraph on this page
   * - ``GetInterfaceChildFaceIds`` (table
       overload) or split tables
     - ``InterfaceChildFaceIdsTest``
     - `Interface split helpers`
   * - ``ValidateGlobalFaceMapping(s)``
     - ``GlobalFaceMappingTest``
     - `Interface split helpers`
   * - ``RegisterLocalInterfaceFace`` / CreateLink mapping
     - ``RegisterLocalInterfaceFaceTest``
     - `Interface split helpers`
   * - ``ExpandBcRecordForSplitInterfaces`` / BC reset on split
     - ``ExpandBcRecordSplitTest``
     - `Interface split helpers`
   * - ``BuildInverseInterfaceMapForZone`` / ConstructNewInterfaceMap
     - ``BuildInverseInterfaceMapTest``
     - `Interface split helpers`
   * - quad diagonal local tris / node map
     - ``QuadDiagonalSplitTest``
     - `Interface split helpers`
   * - ``GridElem`` scan/export/reorder
       orchestration only
     - tests for any moved pure step; full
       CGNS path only if behavior changes
     - `Pipeline overview` if steps change

Suggested PR description lines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy and tick what applies::

    [ ] Staged-face or BC write path covered by FaceSolver* tests
    [ ] Reorder / boundary layout covered by FaceTopoReorderLinkTest
    [ ] Identity/orientation or interface helpers covered by named tests
    [ ] doc/source/face_stage_pipeline.rst updated (or explicitly N/A)
    [ ] Open design pressures updated if testability did not improve

Related pages
^^^^^^^^^^^^^

* :doc:`grid_topology` -- runtime topology contracts and derived connectivity
* :doc:`mesh_architecture` -- mesh information levels and adapters

Source anchors
--------------

* ``codes/geometry/include/FaceSolver.h``
* ``codes/geometry/include/FaceTopo.h``
* ``codes/geometry/include/FaceSearch.h``
* ``codes/geometry/include/IFaceLink.h``
* ``codes/geometry/src/GridElem.cpp``
