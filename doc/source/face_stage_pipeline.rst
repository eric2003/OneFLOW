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

Polyhedral-face staging and deduplication
-----------------------------------------

``GridElem::ScanPolyhedralFaces`` forwards only ``NGON_n`` sections to
``FaceSolver::StagePolyhedralFacesFromConnectivity``. The staging algorithm has a
validation phase followed by a mutation phase:

#. Validate the offset table has one boundary offset per polyhedral face plus one.
#. Validate every face span is ordered and lies within the connectivity array.
#. Reject any 3D polyhedral face with fewer than three node references.
#. Only after all spans pass, build each node ring and query ``faceLookup``.
#. Call ``StageNewFace`` only for a new undirected identity; duplicate rings
   reuse the existing staged face.

Validating the full section before lookup or topology mutation gives malformed
later faces a strong no-partial-staging guarantee. Deduplication identifies a
face by its node set; the stored ring remains the first-seen ordering used by
later orientation/ownership logic. The current staging contract checks the
number of node references for a 3D polyhedral face but does not reject repeated
node IDs within a ring. This minimum does not apply to 2D edge faces.

This routine is valid only for 3D polyhedral face connectivity. The current
``GridElem::HasPolygonSection`` dispatch still treats either ``NGON_n`` or
``NFACE_n`` as a signal to enter this path; dimension-aware classification of
NGON-only sections is a documented adapter gap in :doc:`grid_topology`.

**Tests:** ``FaceSolverStageTest.StagePolyhedralFacesRejectsFewerThanThreeNodes``
and ``FaceSolverStageTest.StagePolyhedralFacesValidatesAllSpansBeforeStaging``.
``FaceSolverStageTest.TwoDimensionalTriangleCreatesTwoNodeEdgeFaces`` locks the
separate 2D cell-to-edge path.
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
  local interface face joins a global id, and rejects duplicate registration
  of the same zone/local-face pair (**test:**
  ``RegisterLocalInterfaceFaceTest``).
* ``FindRemoteInterfacePartner`` selects a non-self zone/local pair only
  when the exact caller (zone and local face) is registered on that global
  face. If the caller is absent, it returns ``false`` without changing output
  arguments, rather than returning an unrelated participant (**test:**
  ``FindRemoteInterfacePartnerTest``). ``MatchInterfaceTopology`` resets
  the stored remote zone and local-face slots to ``-1`` before each lookup, so
  a rematch cannot retain a partner from an earlier mapping state.
  ``MatchPeriodicInterface`` preserves the existing topological partner for
  faces with exactly two zone participants. For other faces, it clears both
  partner fields before geometric matching; a failed point/face lookup therefore
  leaves the face explicitly unmatched instead of retaining a stale partner.\n  The periodic partner row is scanned for the first participant other than the\n  current zone/local-face pair. A self-only row leaves the face unmatched, and\n  mismatched zone/local row lengths are rejected (**tests:**\n  ``FindPeriodicInterfacePartnerTest``).
* ``CountInterfaceFacesMissingPartner`` counts local interface faces whose
  global face has fewer than two zone slots (**test:**
  ``CountInterfaceFacesMissingPartnerTest``). The runtime diagnostic is
  therefore labeled "Interface faces missing a partner", not "Periodic boundary
  faces": this counter examines the interface mapping tables and does not
  independently classify periodic boundaries. Faces with more than two zone
  slots remain excluded for legacy compatibility; that policy is not changed
  by the diagnostic wording correction.
* ``ExpandBcRecordForSplitInterfaces`` repeats BC slots when an interface face
  splits into children (**test:** ``ExpandBcRecordSplitTest``).
* ``BuildInverseInterfaceMapForZone`` fills ``gI2ZidNew`` / ``g2lNew`` from a
  zone's ``l2gNew`` row (at most two zone slots). It validates all global-face
  IDs and existing referenced rows before mutating either inverse table, so an
  invalid later ID cannot leave a partial update (**tests:**
  ``BuildInverseInterfaceMapTest``).
* ``RebuildZoneL2gNewAfterInterfaceSplit`` rebuilds one zone's ``l2gNew`` /
  ``nChild`` when parents split (**test:** ``RebuildZoneL2gNewTest``).
* ``RebuildFaceCellSidesAfterInterfaceSplit`` expands left cells and assigns
  ghost right cells after interface splits (**test:**
  ``RebuildFaceCellSidesTest``).
* ``CollectBoundaryInterfaceChildFaceIds`` walks boundary BC types and records
  per-interface child global-id lists; ``CollectBoundaryInterfaceChildCounts``
  is the size-only view (**test:**
  ``CollectBoundaryInterfaceChildCountsTest``).
* ``RebuildFaceNodeListsAfterInterfaceSplit`` rebuilds face-to-node lists after
  interface splits (**test:** ``RebuildFaceNodeListsTest``).
* Quad diagonal candidates:
  ``BuildQuadDiagonalLocalTris`` / ``MapLocalTrisToGlobalNodes``
  (**test:** ``QuadDiagonalSplitTest``).


Interface mapping names and design rules
-----------------------------------------

The current names are historical abbreviations, not self-explanatory concepts.
Use the following interpretation when reading or changing this code:

.. list-table:: Interface mapping vocabulary
   :header-rows: 1
   :widths: 24 36 40

   * - Existing name
     - Concept
     - Index meaning
   * - ``l2g``
     - local-to-global map
     - ``l2g[zone][localFace]`` is a global interface-face ID
   * - ``gI2Zid``
     - global Interface-face to zone IDs (the ``I`` means Interface)
     - ``gI2Zid[globalFace][slot]`` is a zone ID
   * - ``g2l``
     - global-to-local map
     - ``g2l[globalFace][slot]`` is a local interface-face ID
   * - ``*New``
     - candidate rebuild state
     - buffers not yet committed to the live mapping
   * - ``nChild``
     - legacy split bookkeeping
     - per-zone counts/slots used by interface-face split rebuilding

The important invariant is more precise than the abbreviations suggest:
``gI2Zid[globalFace][slot]`` and ``g2l[globalFace][slot]`` form one
paired reference, so their row lengths and slot order must agree. A local
entry in ``l2g`` points back to the global face, but the current validator
checks the paired inverse tables' shape; it does not prove every semantic
round-trip across all three tables.

Naming is a design tool, not a contest to make every identifier longer.
Prefer names that communicate the *role at the current abstraction level*:
`globalFaceId`, `zoneId`, and `localFaceId` are useful scalar names, while
the table names should describe their mapping direction or relation. Avoid
mechanical expansions that repeat context already clear from the type or
function name.

Do not rename these public ``IFaceLink`` fields in one broad sweep. They are
used across topology construction and rebuild paths, and a mass rename would
create churn without simplifying the algorithm. The incremental migration
policy is:

#. Use descriptive local names in new code and pure helpers.
#. Keep the mapping contract documented beside the declarations and tests.
#. When a stable owner for the mapping state is introduced, add semantic
   accessors or a small mapping type and migrate call sites in focused steps.
#. Remove legacy field names only after all production and test call sites have
   moved and the public compatibility impact is understood.

This policy is part of the mesh development handbook: a rename is justified
when it reduces the mental model or prevents a real class of mistakes, not
merely because an abbreviation looks old.

**Current capacity caveat:** ``BuildInverseInterfaceMapForZone`` stores at most
two zone/local references for a global interface face. The current implementation
silently ignores a third reference, and
``BuildInverseInterfaceMapTest.DoesNotExceedTwoZoneSlots`` records that
behavior. This is not equivalent to validating that the input topology is valid.
The two-slot rule appears to model one interface pair; before changing the
behavior to reject a third participant, confirm whether multi-zone junctions
are excluded by the mesh model and by real input cases. If they are excluded,
a clear error may be safer than silent truncation. Until that contract is
confirmed, preserve the current behavior and treat it as an explicit design
question rather than assuming the test proves the policy is correct.

Interface mapping algorithm
---------------------------

Multi-zone meshes share interface faces. OneFLOW keeps **two dual views** of
the same relation:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Table
     - Meaning
   * - ``l2g[zone][localIface]``
     - Local interface face index to global interface face id
   * - ``gI2Zid[global]`` / ``g2l[global]``
     - Global face to the list of (zone, local) participants (usually two)

**Registration (scan time).** Each zone registers its interface faces with
``RegisterLocalInterfaceFace`` after a geometry identity lookup
(``faceLookup`` / ``FaceSearchBasic``). A new global id appends parallel rows
to ``gI2Zid`` and ``g2l``; a shared face pushes a second (zone, local) pair.

**Partner selection.** After all zones register, ``FindRemoteInterfacePartner``
picks the (zone, local) entry on a global face that is not the caller. Faces
with fewer than two zone slots are counted by
``CountInterfaceFacesMissingPartner`` (legacy: more than two slots are not
counted as "missing").

**Split rebuild.** Interface parents may split (for example quad to triangles).
Rebuild runs on ``New`` buffers:

#. ``InitNewInterfaceMappingBuffers`` -- copy inverse maps, empty ``l2gNew``
#. Per zone: ``RebuildZoneL2gNewAfterInterfaceSplit`` -- parents become children
   in ``l2gNew`` / ``nChild``
#. ``ExpandBcRecordForSplitInterfaces`` / face-node and cell-side rebuilds
#. ``BuildInverseInterfaceMapForZone`` -- validate the complete zone row and
   existing inverse rows, then fill ``gI2ZidNew`` / ``g2lNew`` (at most two zone
   slots per global face)
#. ``CommitInterfaceMappingBuffers`` -- validate then swap ``New`` into live tables

The dual-table invariant after commit is the same as after registration:
``ValidateGlobalFaceMappings`` requires equal top-level sizes and equal per-face
list lengths for ``gI2Zid`` and ``g2l``.

Commit is transactional with respect to the live tables: validation runs
against the ``New`` inverse maps before any swap. If validation throws, the
current ``gI2Zid``, ``g2l``, and ``l2g`` remain unchanged and the candidate
buffers remain available for diagnosis or correction. The regression test
``CommitInterfaceMappingBuffersTest.FailedValidationLeavesLiveTablesUnchanged``
locks this failure behavior.


Boundary-condition assignment algorithm
---------------------------------------

Staged faces start with BC slots set to ``INVALID_INDEX``. CGNS BOCO regions
supply a **vertex set** (nodes on the boundary patch). Assignment walks staged
faces and writes BC only through ``SetStagedFaceBc``.

**Match predicate.** ``FaceNodesContainedInVertexSet`` (used by
``CheckBcFace``) is true when every node of a face appears in the BOCO vertex
set. The BOCO set may contain extra nodes; face node order does not matter.

**Eligibility filter.** Before the vertex-set match,
``IsEligibleForStagedBcAssignment`` decides whether a face may be written:

#. ``currentBcType`` must still be ``INVALID_INDEX`` (first writer wins; later
   BOCO regions do not overwrite).
#. If ``boundaryFacesOnly`` is true, ``rightCell`` must be ``INVALID_INDEX``
   (open boundary). Detail mode drops this restriction.

**Two scan modes.**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - API
     - Filter
   * - ``ScanBcFace``
     - ``boundaryFacesOnly=true``: open right cell and no BC yet
   * - ``ScanBcFaceDetail``
     - no BC yet only; may label interior faces if nodes lie in the set

Both call ``AssignBcToMatchingFaces`` then ``FaceNodesContainedInVertexSet``.

**Interface fill.** ``ScanInterfaceBc`` assigns ``BCTypeNull`` (and name id
``-1``) to any remaining face still at ``INVALID_INDEX``. That closes the
staged BC vector before export.

**Export.** ``CollectBoundaryBcEntries`` (used by
``AppendBoundaryBcToRecord``) walks faces in **current order** and appends
staged BC type/name for every face with ``right == INVALID_INDEX``. Interior
faces are skipped. The helper **appends** (does not clear) the output arrays,
matching historical ``BcRecord`` accumulation.

After ``ReorderLink``, boundary faces form a prefix and relative boundary
order is preserved, so ``BcRecord`` stays aligned without a second
permutation (``ValidateBoundaryBcLayout``).

**Tests:** ``FaceNodesContainedInVertexSetTest``,
``CollectBoundaryBcEntriesTest``, ``FaceSolverStageTest``,
``FaceTopoReorderLinkTest``.


Face ownership algorithm
------------------------

Each staged face has left/right cell slots and a claim flag used by the polygon
path.

**First touch (new face).** ``StageNewFace`` records the creating cell as left
owner when provided; the face flag starts as ``kFaceUnclaimed`` for polygon
staging, or is claimed when volume elements stage with a known left cell.

**Polygon claim path.**

* ``ClaimStagedFaceAsLeft`` -- requires ``kFaceUnclaimed``; sets flag to
  ``kFaceClaimed``, left = cell, right = ``INVALID_INDEX``, clears staged BC.
* ``AttachStagedFaceAsRight`` -- requires already claimed; sets right = cell.

**Volume element path.** ``ScanElementFace`` uses geometry identity
(``faceLookup``). New faces go through ``StageNewFace``. Existing faces call
``BindCellToExistingFace``, which applies ``BindCellToFaceSides``:

#. If left is open (``INVALID_INDEX``), assign left = cell.
#. Else if right is open and left is a different cell, assign right = cell.
#. Else no change (already bound, or the same cell is already on left).

This implements the usual unstructured-mesh rule: the first incident cell owns
the left side; a second distinct cell becomes the right (interior) neighbour.
A second visit by the same cell does not fill the right slot (avoids marking a
boundary face as interior by mistake when orientation/duplication appears).

**Tests:** ``BindCellToFaceSidesTest``, ``FaceSolverStageTest``.


Staged parallel-array growth
----------------------------

``StageNewFace`` is the only growth path for a newly discovered unique face.
One call appends **one** entry to each of:

* topology: ``faces``, ``faceTypes``, ``faceFlags``
* ownership: ``leftCells``, ``rightCells`` (right starts as ``INVALID_INDEX``)
* local indices: ``lPosition``, ``rPosition`` (seeded to 0)
* BC: ``faceBcTypes``, ``faceBcNameIds`` (both ``INVALID_INDEX``)
* split: ``childFid`` (empty list)

After growth, ``StagedBcMatchesFaceCount`` requires all of these lengths to be
equal. That invariant is what makes later pure helpers safe: BC export,
ownership bind, and child counting never re-size topology.

**Residual interface BC.** After BOCO vertex matching, any face still at
``INVALID_INDEX`` is treated as an interface (or unset) face.
``FillUnsetStagedBcSlots`` (used by ``ScanInterfaceBc``) writes
``BCTypeNull`` and name id ``-1`` into those slots only. Already assigned BC
types are left unchanged.

**Simple faces.** ``CountEmptyChildFaceLists`` counts ``childFid`` rows with
no children; ``CountFacesWithoutChildren`` checks alignment with ``faces``
then delegates.

**Tests:** ``FillUnsetStagedBcSlotsTest``, ``FaceSolverStageTest``.


Boundary-first reorder algorithm
--------------------------------

After staging and BC export, faces may still be in discovery order (volume
sweep). Solvers expect **all boundary faces in a contiguous prefix**.

``BuildBoundaryFirstPermutation(rightCells)`` builds ``newToOld``:

#. First pass: append indices with ``right == INVALID_INDEX`` (boundary).
#. Second pass: append the remaining indices (interior), preserving relative
   order within each class (stable partition).

``ReorderLink`` applies that permutation to every face-indexed array
(faces, types, flags, left/right cells, positions, and related metadata) and
clears derived ``c2f``. ``BcRecord`` is **not** reordered separately: it was
filled in discovery order of boundary faces, and the partition preserves the
relative order of those boundary faces, so the prefix of reordered faces still
matches ``BcRecord`` (``ValidateBoundaryBcLayout``).

**Tests:** ``FaceTopoReorderLinkTest``.

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
   * - ``StagePolyhedralFacesFromConnectivity``
     - ``FaceSolverStageTest.StagePolyhedralFacesRejectsFewerThanThreeNodes``,
       ``FaceSolverStageTest.StagePolyhedralFacesValidatesAllSpansBeforeStaging``
     - `Polyhedral-face staging and deduplication`
   * - ``FillUnsetStagedBcSlots`` / ScanInterfaceBc /
       ``CountEmptyChildFaceLists``
     - ``FillUnsetStagedBcSlotsTest``
     - `Staged parallel-array growth`
   * - ``SetStagedFaceBc`` /
       ``AssignBcToMatchingFaces`` /
       BC slot privacy or accessors
     - ``FaceSolverStageTest`` (BC write/read)
     - write-path table on this page
   * - ``FaceNodesContainedInVertexSet`` / CheckBcFace
     - ``FaceNodesContainedInVertexSetTest``
     - `Boundary-condition assignment algorithm`
   * - ``IsEligibleForStagedBcAssignment``
     - ``IsEligibleForStagedBcAssignmentTest``
     - `Boundary-condition assignment algorithm`
   * - ``CollectBoundaryBcEntries`` / AppendBoundaryBcToRecord
     - ``CollectBoundaryBcEntriesTest``
     - `Boundary-condition assignment algorithm`
   * - ``ClaimStagedFaceAsLeft`` /
       ``AttachStagedFaceAsRight`` /
       ``BindCellToExistingFace``
     - ``FaceSolverStageTest`` (ownership)
     - post-stage writes table on this page
   * - ``BindCellToFaceSides``
     - ``BindCellToFaceSidesTest``
     - `Face ownership algorithm`
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
   * - ``FindRemoteInterfacePartner`` / MatchInterfaceTopology
     - ``FindRemoteInterfacePartnerTest``
     - `Interface split helpers`
   * - ``CountInterfaceFacesMissingPartner``
     - ``CountInterfaceFacesMissingPartnerTest``
     - `Interface split helpers`
   * - ``InitNewInterfaceMappingBuffers`` / Commit / InitNewLgMapping
     - ``InitNewInterfaceMappingBuffersTest``
     - `Interface mapping algorithm`
   * - ``ExpandBcRecordForSplitInterfaces`` / BC reset on split
     - ``ExpandBcRecordSplitTest``
     - `Interface split helpers`
   * - ``BuildInverseInterfaceMapForZone`` / ConstructNewInterfaceMap
     - ``BuildInverseInterfaceMapTest``
     - `Interface split helpers`
   * - ``RebuildZoneL2gNewAfterInterfaceSplit`` / ModifyBoundaryInformation
     - ``RebuildZoneL2gNewTest``
     - `Interface split helpers`
   * - ``RebuildFaceCellSidesAfterInterfaceSplit`` / SetNewFace2Cell
     - ``RebuildFaceCellSidesTest``
     - `Interface split helpers`
   * - ``CollectBoundaryInterfaceChildFaceIds`` / Counts
     - ``CollectBoundaryInterfaceChildCountsTest``
     - `Interface split helpers`
   * - ``RebuildFaceNodeListsAfterInterfaceSplit`` / SetNewFace2Node
     - ``RebuildFaceNodeListsTest``
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
