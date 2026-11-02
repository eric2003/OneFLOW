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


Why this local refactor exists
------------------------------

Architecture map (local refactor status)
----------------------------------------

This is the intentional end state of the **current** local refactor, not a
wishlist for the whole solver.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Area
     - Closed
     - Still open
   * - Face discovery / stage
     - ``StageNewFace``, Claim/Attach/Bind,
       BC assign/export, boundary-first reorder
     - End-to-end CGNS zone fixtures
   * - Interface dual maps
     - Register / Init / Commit pure helpers;
       zone bookkeeping accessors;
       dual-map fields private;
       multi-zone split e2e composition
       (``MultiZoneInterfaceSplitE2ETest``)
     - Whether to strip orphan inverse rows
       (detection: ``CollectUnreferencedGlobalInterfaceFaces``)
   * - IFaceLink services
     - ``faceLookup``, ``point_search``,
       ``face_search`` private + accessors;
       injectable ``PointLocator``
       (``ReplacePointLocator``)
     - Production ADT still needs nodeMesh
       for non-injected paths
   * - Child mapping after split
     - Parent-specific ``ByParent`` tables;
       two-parent shared-child test;
       multi-zone split e2e without CGNS;
       pure line midpoint + candidate maps
     - Production still needs locator +
       zone ``nodeMesh`` (pure prefix locked by
       ``LineSplitPureStagesTest``)

**Rule of thumb:** if a change only renames a pure predicate already covered
by tests, skip it. Prefer closing an open row in this table.

What this local refactor is for
-------------------------------

The geometry stack once mixed **discovery**, **ownership**, and **physics labels**
in the same loops, with multi-zone interface tables as public data bags. The
local refactor is not a rewrite of the CFD solver; it is a **boundary-drawing**
exercise on the face/interface path so that:

#. Each job has a named write path (stage/claim/bind, dual-map Register/Init/Commit,
   BC export + boundary-first reorder as one transaction).
#. Invariants are checkable in ``tests/geometry`` without a full CGNS case.
#. Hard edges (ADT point search, multi-zone split rebuild) become ports with
   accessors or injection, not more public fields.

If a change does not make a contract easier to state or test, it is probably
detail noise. Prefer closing a row in the architecture map over renaming a
stable pure predicate.



The face / interface / BC path used to interleave **three different jobs** in
the same loops:

#. **Discover** unique faces while scanning elements (geometry identity).
#. **Own** each face from the left and right cells (topology).
#. **Label** boundary faces from CGNS BOCO vertex sets (physics BC).

Those jobs share parallel arrays (``faces``, left/right cells, BC slots,
``l2g`` / ``gI2Zid``). Before the staged APIs, alignment was accidental: each
scan path did its own ``push_back``. Unit tests had to stand up a full CGNS
zone to assert a single rule.

The local refactor's meaning is therefore not "more free functions". It is:

* **One growth contract** -- ``StageNewFace`` and
  ``StagedBcMatchesFaceCount`` so parallel arrays stay aligned by construction.
* **Named write paths** -- Claim / Attach / Bind / ``SetStagedFaceBc`` so
  ownership and BC are not ad-hoc index writes.
* **Pure rules extracted from adapters** -- eligibility, vertex-set match,
  boundary-first permutation, interface dual maps -- so intent is unit-tested
  without CGNS.
* **Documented algorithms** on this page so the next change can see the
  *pipeline*, not only the last helper.

What remains hard (see `Open design pressures`_) is still the right backlog:
``IFaceLink`` public tables, ``PointLocator`` injection, parent-relative child
mappings after interface split, and true end-to-end CGNS fixtures. Prefer
closing those seams over extracting another one-line predicate.

GridElem orchestration
----------------------

``GridElem`` is the **conductor** for one zone's calculation grid. It does not
implement face identity or BC matching; it sequences already-tested steps:

#. Scan elements into ``FaceSolver`` (stage / claim / bind).
#. Scan BOCO into staged BC (``AssignBcToMatchingFaces`` / residual fill).
#. ``AppendBoundaryBcToRecord`` -- export boundary BC in discovery order.
#. Move ``FaceTopo`` onto ``UnsGrid``, map CGNS BC ids
   (``BcTypeMap::MapFromCgnsInPlace``).
#. ``ReorderBoundaryFirstAndValidate`` -- boundary-first face order and
   ``BcRecord`` layout check as **one** post-stage transaction.

If a rule can be stated without ``GridElem`` or ``CgnsSection``, it belongs in a
pure helper and a geometry unit test. ``GridElem`` should stay thin.

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

The split table uses two index domains because face identity and
parent-relative reconstruction are different concepts:

* ``childFaceIdsOfParent[globalParentFaceId]`` lists global child face IDs.
* ``childFaceNodeIds[globalChildFaceId]`` and
  ``childFaceNodeFlags[globalChildFaceId]`` are legacy child-indexed
  mapping rows retained for compatibility.

Candidate resolution and local mapping
--------------------------------------

For a quad, ``BuildQuadDiagonalLocalTris`` emits four candidate triangles in
two possible diagonal decompositions. These candidates are ordered by the
algorithm, not by which faces happen to have been registered globally.
``ResolveRegisteredChildFaces`` therefore returns a pair for each match:
``globalFaceId`` identifies the registered face, while ``candidateIndex``
points back to the original candidate triangle.

The distinction matters for partial matches. If only candidate 3 is registered,
it is still candidate 3 for local-node and flag lookup; it must not become
candidate 0 merely because it is the first item in the filtered result. The
C++ integration tests ``QuadSplitUsesOriginalCandidateIndexForPartialMatch``
and ``QuadSplitPreservesLastCandidateIndexForPartialMatch`` check this through
``FaceSearch::SplitQuad2Tri`` and its published parent-relative mapping, rather
than testing the resolver in isolation.

This is a small but useful architectural boundary: identity resolution decides
which global children exist; the original candidate index preserves the
relationship between that identity and the geometric decomposition. Keeping
both values together prevents a filtered list from silently changing the
meaning of a local mapping row.
* ``childFaceNodeIdsByParent[globalParentFaceId]`` and
  ``childFaceNodeFlagsByParent[globalParentFaceId]`` store the preferred
  mappings aligned by child position within that parent's child list.

``ResetForFaceCount`` sizes the global-ID-indexed tables to the global face
count. Parent-specific rows are aligned to ``childFaceIdsOfParent[parent]``,
not indexed by child ID. A flag of ``1`` means the corresponding node ID is
a corner index in the parent face ring; other flag values retain the mapping
semantics defined by the split builder.

This distinction is required when a global child face is shared by multiple
parents. The child ID can be shared because face identity is order-independent,
but local corner indices belong to a particular parent ring. For example, the
same triangle with global nodes `B, C, D` can be described by local corners
`[1, 2, 3]` in a parent ring `[A, B, C, D]` and `[0, 1, 2]` in a
different parent ring `[B, C, D, E]`. Reusing one child's
mapping row for both parents would reconstruct the wrong ordered node ring for
at least one consumer.

The production path therefore prefers parent-specific rows when collecting
node mappings for interface rebuild. Legacy global-child rows are a fallback
only when both parent-specific rows are absent or empty. If either new row is
present but the pair does not align with the parent's child count, collection
throws instead of silently masking malformed new-format metadata. An empty
child list is not an error: it means the parent was not split, so the rebuild
preserves that face without requiring node mapping metadata.

* ``GetInterfaceChildFaceIds(l2g, splitTable, zone, local)`` resolves
  parent-to-child identity and validates the local and parent lookup indices.
  It cannot validate child IDs against the global face count because the split
  table does not own that count. It deliberately does not inspect node mapping
  rows; ``CollectBoundaryInterfaceChildNodeMappings`` selects parent-specific
  rows or the legacy fallback and validates the selected representation. It
  stages both output tables locally and publishes them only after every
  interface entry is valid, so a later malformed parent cannot expose a
  partial result (**tests:** ``InterfaceChildFaceIdsTest``,
  ``CollectBoundaryInterfaceChildCountsTest``, and
  ``RebuildFaceNodeListsTest.MappingCollectorDoesNotPublishPartialOutputs``).
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
  leaves the face explicitly unmatched instead of retaining a stale partner.
  The periodic partner row is scanned for the first participant other than the
  current zone/local-face pair. A self-only row leaves the face unmatched, and
  mismatched zone/local row lengths are rejected (**tests:**
  ``FindPeriodicInterfacePartnerTest``).
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
* ``RegisterLocalInterfaceFace`` constructs the initial global-face participant
  rows. A global face supports at most two zone participants; attempting to
  register a third participant now throws before changing ``l2g``, ``gI2Zid``,
  or ``g2l`` (**test:** ``RegisterLocalInterfaceFaceTest``).
* ``BuildInverseInterfaceMapForZone`` refreshes the inverse references for
  one zone from its rebuilt ``l2gNew`` row. Since ``InitNewInterfaceMappingBuffers``
  starts from copies of the live inverse maps, a rebuild must update the
  existing zone slot's local-face index rather than append that zone a second
  time. If the zone is absent, it is appended only when the global face has
  fewer than two participants. The legacy policy still silently ignores a
  third distinct participant rather than expanding the data model. This differs
  from initial registration, which rejects unsupported topology early. The
  routine validates all global-face IDs, per-zone uniqueness, and existing
  referenced rows before mutation (**tests:**
  ``BuildInverseInterfaceMapTest.RefreshesExistingZoneSlotDuringRebuild``,
  ``BuildInverseInterfaceMapTest.DoesNotDuplicateExistingZoneParticipant``).
* ``RebuildZoneL2gNewAfterInterfaceSplit`` rebuilds one zone's ``l2gNew`` /
  ``nChild`` when parents split. It validates all parent and child global IDs
  before clearing either output, so a negative ID cannot leave a partial
  candidate map (**tests:** ``RebuildZoneL2gNewTest``).
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
  ``BuildQuadDiagonalLocalTris`` enumerates four candidates, not one selected
  triangulation: rows 0-1 use diagonal 0-2 and rows 2-3 use diagonal 1-3.
  For local corners ``0, 1, 2, 3``, the expected pairs are
  ``(0,1,2) + (2,3,0)`` and ``(3,0,1) + (1,2,3)``. Each pair must
  preserve the four boundary edges and share its diagonal exactly once per
  triangle. ``FaceSearch::SplitQuad2Tri`` then maps all four candidates to
  global node IDs and keeps only triangles already found by ``FindFace`` in
  the global face registry. Thus the runtime result depends on the registry:
  when only one diagonal is registered, exactly its two triangle faces should
  be recorded (test: ``FaceSearchBasicTest.QuadSplitUsesOnlyRegisteredTrianglesForOneDiagonal``).
  ``FindFace`` identifies a face by its node set, independent of node order,
  while ``GetFaceNodeRing`` preserves the first registered ordering. The runtime
  split must therefore find an existing child even when its stored ring differs
  from the candidate order (test:
  ``FaceSearchBasicTest.QuadSplitMatchesRegisteredTrianglesRegardlessOfNodeOrder``).
  Zero or one points to missing/inconsistent face registration, while more than
  two means both candidate diagonals may be represented and needs topology-level
  investigation. The candidate builder itself does not choose a diagonal or
  compare triangle quality.

  Diagnostic order for a failed quad split:
  
  #. Check the parent ring has exactly four nodes and inspect its stored order.
  #. Check the four local candidate rows and confirm each pair shares one
     diagonal while retaining the boundary edges.
  #. Check local-to-global mapping against the parent ring.
  #. Query the face registry for each mapped candidate; distinguish a missing
     triangle from a mapping mismatch.
  #. Only after these checks pass, inspect child registration and the later
     interface-map / boundary / cell-side rebuilds.

  **Tests:** ``QuadDiagonalSplitTest.BuildsFourCandidateTrisForBothDiagonals``,
  ``QuadDiagonalSplitTest.EachDiagonalCandidatePairPreservesBoundaryAndSharesOneDiagonal``,
  and ``QuadDiagonalSplitTest.MapsLocalCornersThroughQuadNodes``.
  ``MapLocalTrisToGlobalNodes`` validates every candidate row has exactly three
  corners and every corner index is in range before replacing the output, so
  malformed later rows cannot leave a partial mapping (**tests:**
  ``QuadDiagonalSplitTest.MapsLocalCornersThroughQuadNodes`` and
  ``FaceSearchBasicTest.InvalidQuadCandidateLeavesMappedOutputUnchanged``).
* ``SplitQuad2Tri`` clears the selected parent's previous child list before
  checking its shape; a direct call on a non-quad face therefore cannot preserve
  a stale split result (**test:** ``FaceSearchBasicTest.NonQuadSplitClearsStaleParentMapping``).
* ``InterfaceFaceSplitTable::ClearParentMapping`` and
  ``SetChildMapping`` reject misaligned parent-specific ID/flag table
  capacities before mutating the parent's child list (**tests:**
  ``FaceSearchBasicTest.ClearParentMappingRejectsMisalignedParentTables``
  and ``FaceSearchBasicTest.SetChildMappingRejectsMissingParentSpecificRow``).


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

Line midpoint split algorithm
-----------------------------

The 2D interface split path separates geometry discovery from mapping assembly.
``FaceSearch::SplitLine`` first clears the selected parent's old child list and
ignores faces that are not two-node edges. ``TrySplitLineAtMidpoint`` validates
the parent face's zone association before querying geometry, then checks whether
the geometric midpoint is already present in ``PointLocator``; this path does
not independently decide to split every edge. Invalid face-to-zone associations
therefore fail before midpoint lookup or zone-local node mutation.

When a midpoint is available, the algorithm resolves both candidate child
edges against the global ``FaceSearch`` registry before touching the zone mesh.
If neither child exists, it returns without appending a zone-local node. The
resolved child IDs and candidate positions are retained and reused when parent
mappings are published, so the registry lookup is performed once per candidate
rather than once for the preflight and again for publication. Duplicate child
identities are rejected by the shared resolver before any zone-local mutation.
``ResolveRegisteredChildFaces`` is now an explicit topology-resolution step
that can be tested without ``IFaceLink`` or a real grid. It skips candidate
edges absent from the registry and returns records that keep each global child
face ID paired with its original candidate position. A single record prevents
the caller from accidentally combining one child's ID with another candidate's
local mapping by indexing two independent output arrays. Returning the records
by value also makes the result all-or-nothing if duplicate identities are
detected. The current policy intentionally permits a partial match: one
registered child is enough to proceed, while no registered children means no
local midpoint insertion. Tests make this policy explicit so a later change to
require both children would be a deliberate behavior change.

The ``registered_child_partial_match`` workbench scenario mirrors
``FaceSearchBasicTest.RegisteredChildResolutionPreservesCandidateIndexForPartialMatch``.
Its input makes the distinction visible: candidate 0 is missing, candidate 1
resolves to global child face 0, and the returned candidate index must remain 1.
The scenario links both the helper-level test and
``FaceSearchBasicTest.QuadSplitUsesOriginalCandidateIndexForPartialMatch``.
The latter runs the real quad-split publication path with only candidate 1
registered and verifies that its parent-local node mapping is ``[2, 3, 0]``.
This catches the original class of bug at the consumer boundary: even if the
resolver reports the right candidate index, the caller must use that index when
publishing parent-relative metadata. The current GUI scenario backend is still
fixture-backed; neither linked C++ test is executed by the GUI. A native test
adapter that passes scenario parameters to C++ remains the next integration
boundary.

The algorithm then resolves the midpoint in two index spaces:

#. ``PointLocator::FindPoint`` returns the existing **global search** midpoint ID.
   If no midpoint is registered there, this edge is not split.
#. The zone's ``nodeMesh`` is queried using the same coordinate. If the point is
   absent there, it is appended and the returned index becomes the **zone-local**
   midpoint ID. The already-discovered global ID is reused when constructing
   child edge identities; the locator is not searched or modified a second time.

The line split now has two deliberately separate representations:

* ``BuildGlobalLineSplitCandidates`` creates directed global child-edge identities
  ``(A, M_global)`` and ``(M_global, B)``. These candidates are resolved
  against the global face registry before the zone mesh is mutated.
* ``BuildLocalLineSplitMapping`` creates zone-local edge indices
  ``(0, M_local)`` and ``(M_local, 1)``, plus endpoint flags
  ``(1, -1)`` and ``(-1, 1)``.

The global midpoint ID identifies topology; the local midpoint ID addresses one
zone's node array. Keeping these outputs separate makes the preflight lookup
independent of the zone-local index and prevents a temporary local-mapping
builder from also producing a global representation that its caller discards. The endpoint flags encode which
original endpoint each child retains. Reversing the parent ring must reverse
the global child-edge direction while leaving the local index pattern and
flags unchanged.

The distinction matters because global face identity is order-independent, but
the stored face ring and local endpoint metadata are directional.
``BuildGlobalLineSplitCandidates`` replaces its output candidate list on each
call, so a reused list cannot retain child edges from an earlier split.

**Tests:** ``FaceSearchBasicTest.GlobalLineCandidatesPreserveEdgeDirection``,
``FaceSearchBasicTest.GlobalLineCandidatesPreserveReversedEdgeDirection``,
``FaceSearchBasicTest.GlobalLineCandidatesReplacePreviousOutput``, and
``IFaceLinkPointLocatorInjectionTest.SplitLineWithoutRegisteredChildrenDoesNotPublishMapping``.

The injected-locator integration test exercises the real ``SplitLine`` entry
point for the negative topology case: a midpoint exists in global search, but
neither child edge is registered. It verifies the algorithm exits before asking
for a zone grid and leaves the parent child list and parent-relative mapping
rows empty. This is more valuable than another helper-only test because it
locks the order of operations across locator injection, global topology lookup,
and mapping publication.

``FaceSearchBasicTest.RegisteredChildResolutionPreservesCandidateIndexForPartialMatch``
locks the complementary table-level rule: when only candidate 1 is registered,
the resolver returns candidate index 1 rather than compacted index 0. That index
selects the correct parent-local node/flag row during publication. The test does
not claim to exercise zone-local midpoint insertion; that remains covered only
by the full ``SplitLine`` path when a suitable grid fixture is available.

**Current boundary:** this regression covers the no-child short-circuit, not a
successful zone-local insertion. The positive path still depends on
``IFaceLink``/``nodeMesh`` services and should be tested only when a small
grid fixture can assert the local node insertion and published mappings
together. Do not add another abstraction just to mock calls that already have
an injection seam.

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
first verifies that the caller's exact (zone, local-face) pair is present in
the global-face row, then selects the first other participant. If the caller
is absent, it returns false without modifying the output values; this prevents
an unrelated participant from being reported as a topological neighbour.

Periodic matching has a different lookup contract. It geometrically transforms
the current face and searches for the corresponding face in the transformed
location. That partner face's row normally belongs to another zone, so the
current (zone, local-face) pair may be absent. ``FindPeriodicInterfacePartner``
therefore selects the first participant unless it is the exact current pair;
a self-only row returns false. Both helpers require aligned zone/local rows and
preserve output values when no partner is found (tests:
``FindRemoteInterfacePartnerTest``, ``FindPeriodicInterfacePartnerTest``).

Faces with fewer than two zone slots are counted by
``CountInterfaceFacesMissingPartner`` (legacy: more than two slots are not
counted as "missing").

**Inverse-map construction.** ``BuildInverseInterfaceMapForZone`` derives
the per-global-face participant rows from one zone's ``l2gNew`` list. A
zone-local list must contain each global face ID at most once: otherwise the
same zone could be inserted twice into one participant row and later be
mistaken for its own remote neighbour. The helper checks all IDs, duplicate
global IDs, and existing row alignment before mutating either inverse table,
so a rejected zone leaves both tables unchanged
(**test:** ``BuildInverseInterfaceMapTest.RejectsDuplicateGlobalFaceWithinZoneWithoutMutation``).

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

``ReorderBoundaryFirstAndValidate`` is the single post-stage entry used by
``GridElem``: reorder then validate. Callers should not reorder without the
layout check.

**Tests:** ``FaceTopoReorderLinkTest``,
``ReorderBoundaryFirstAndValidateTest``.


IFaceLink dual-map boundary
---------------------------

``IFaceLink`` holds the multi-zone interface relation in two directions
(``l2g`` and ``gI2Zid`` / ``g2l``). That is intentional: registration,
partner match, and split rebuild all need both views.

The local refactor does **not** pretend these tables are gone. It draws a
boundary around *how* they change:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Operation
     - Allowed path
   * - Register a local face on a global id
     - ``CreateLink`` -> ``RegisterLocalInterfaceFaceOnMaps``
   * - Stage / commit after interface split
     - ``InitNewInterfaceMappingBuffers`` /
       zone rebuild helpers /
       ``CommitInterfaceMappingBuffers``
   * - Find remote partner / count missing
     - ``FindRemoteInterfacePartner`` /
       ``CountInterfaceFacesMissingPartner``
   * - Geometry identity for interface nodes
     - private ``faceLookup`` inside ``CreateLink`` only

``faceLookup``, ``point_search``, and ``face_search`` are private. Use
``CreateLink``, ``GetPointLocator`` / ``EnsurePointLocatorReady``,
``GetFaceSearch``, and ``MutableInterfaceSplit`` / ``GetInterfaceSplit``.
Unit tests may call ``ReplacePointLocator`` with a ``PointLocator`` subclass
and ``markReady=true`` to avoid building the production ADT from grids.

Split rebuild zone rows go through ``EnsureSplitBookkeepingMatchesZones``,
``GetL2gForZone``, ``MutableL2gNewForZone``, and ``MutableNChildForZone``.
Staging inverse maps use ``MutableGlobalToZoneIdsNew`` /
``MutableGlobalToLocalIdsNew``. ``FaceTopo`` and ``UnsGrid`` read live
``l2g`` through ``GetLocalToGlobalMap`` / ``GetL2gForZone``. Dual-map
fields (``l2g``, ``gI2Zid``, ``g2l``, ``*New``, ``nChild``) are private.

Tests and topology-only clients seed maps with
``EnsureInterfaceZoneCount``, ``EnsureLocalInterfaceFaceCount``, and
``RegisterLocalInterfaceFace`` (member wrapper). Read-only observers use
``GetLocalToGlobalMap``, ``GetGlobalToZoneIds``, and ``GetGlobalToLocalIds``.

**Tests:** ``IFaceLinkMapAccessorsTest``, ``IFaceLinkSplitBookkeepingTest``.


Parent-relative child mapping algorithm
---------------------------------------

Interface split may turn one parent face into several child faces. Child
geometry is stored in two layers:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Table
     - Role
   * - ``childFaceIdsOfParent[parent]``
     - Ordered global child ids for that parent
   * - ``childFaceNodeIdsByParent[parent][i]`` /
       ``childFaceNodeFlagsByParent[parent][i]``
     - Parent-local node encoding for child ``i`` (flag ``1`` = index into
       **this** parent's node ring)
   * - ``childFaceNodeIds[child]`` /
       ``childFaceNodeFlags[child]``
     - Optional global seed (first writer); **not** authoritative when
       parent-specific rows exist

**Why parent-specific rows.** The same undirected child ring can appear under
two parents with different node sets, for example parent
``[A,B,C,D]`` and parent ``[B,C,D,E]`` sharing child ``[B,C,D]``. Local corner
indices differ (``[1,2,3]`` vs ``[0,1,2]``). A single global row keyed only by
child id cannot hold both encodings. Rebuild must read the
**current parent's** mapping row, then ``ResolveChildFaceNodeIds`` against that
parent's ring.

**Publish path.** ``ResolveRegisteredChildFaces`` decides which candidates are
already registered global faces and rejects duplicate child ids for one
parent. ``SetChildMapping`` then writes the parent-specific row (and seeds the
global row only if still empty).

**Consume path.** ``CollectBoundaryInterfaceChildNodeMappings`` uses
parent-specific rows once any such mapping has been published. A split parent
with a missing parent-specific row is an error; falling back to a global child
row could silently reuse another parent's local node numbering for a shared
child. The global child tables remain a compatibility path only when all
parent-specific rows are empty. ``RebuildFaceNodeListsAfterInterfaceSplit``
applies ``ResolveChildFaceNodeIds`` per child.

**Regression:** ``CollectBoundaryInterfaceChildNodeMappingsTest.RejectsMissingParentRowInsteadOfUsingSharedChildGlobalRow``
proves the collector rejects an incomplete modern mapping and leaves output
buffers unchanged. Legacy fixtures that only populate global child rows remain
supported.

**Tests:** ``TwoParentSharedChildMappingTest``,
``FaceSearchBasicTest`` (duplicate / conflict guards).


Geometry workbench (parallel to refactor)
-----------------------------------------

Unit tests lock contracts in CI. A **workbench** is a different surface:

* Operator edits discrete tables (``l2g``, parents, child lists, BC rows).
* Host calls the **same** C++ pipeline the tests call (today:
  ``RunSymmetricInterfaceSplitMappingPipeline``).
* Results are shown as tables (and later drawings); human or golden compare
  judges pass/fail.

This can advance **in parallel** with the geometry refactor:

#. **Fake-first UI:** panels bind to in-memory tables and a stub that returns
   hand-written outputs, so interaction and visualization ship before every
   production edge is free of CGNS.
#. **Contract-ready core:** pure helpers and table pipelines land in
   ``codes/geometry`` (no second algorithm in the UI process).
#. **Swap:** replace the stub with the real pipeline when the scenario matches;
   do not re-implement split logic in Qt/Python.

What the workbench is *not*: a button that only launches GoogleTest. What it
*is*: an editable, visual driver of table-level geometry contracts -- the same
contracts the local refactor is carving out.

Local refactor scorecard
------------------------

Closed enough to stop polishing the same layer:

* Face stage write path and BC export / boundary-first reorder.
* Dual-map mutation boundary (private fields, Register / Init / Commit).
* PointLocator as a replaceable port for tests.
* Multi-zone split **composition** without CGNS.

Still product work, not rename work:

* Production mid-point split still needs real ``nodeMesh`` + ADT unless injected.
* Orphan parent global rows after InitNew copy (document, do not silent-fix).
* Full CGNS zone fixtures remain heavy; prefer more pure compositions first.



Investment priority (do not spread effort evenly)
-------------------------------------------------

Highest return first. Skip lower rows until a higher row is blocked.

#. **Numerical / solver contracts already in motion** (limiter ranges, Venkat
   policy, Marray index/fill). These sit on the hot path of every residual;
   geometry pure-helpers no longer move residual quality.
#. **Geometry hard edges only** — not more one-line predicates:

   * Line/quad split still needs live ``PointLocator`` + zone ``nodeMesh``.
   * Prefer one fixture that exercises that seam (injected locator + minimal
     mesh) over extracting another arithmetic helper.
   * Orphan inverse rows: detect (done); strip only with an explicit product
     decision and migration test.
#. **Workbench → native table pipeline** when development velocity needs it
   (stub is enough until someone edits scenarios daily).
#. **Full CGNS zone e2e** last among geometry items: high confidence, high cost;
   pure compositions and ``LineSplitPureStagesTest`` already cover most
   dual-map and split-table contracts.

**Geometry plateau rule:** if a change only renames or re-wraps a helper already
locked in ``tests/geometry``, stop. The face/interface **boundary-drawing**
goal is largely met; further investment belongs on hard edges or on solver
contracts (priority 1).

``LineSplitPureStagesTest`` locks the pure prefix of line split (sole-zone →
midpoint → candidates → resolve → local mapping) without ADT or CGNS, so the
remaining open map row stays honest: production still needs locator + nodeMesh.


Multi-zone interface split (composed scenario)
----------------------------------------------

Isolated helpers lock one zone or one table. The composed entry point is
``RunSymmetricInterfaceSplitMappingPipeline`` (table in, table out; no
``IFaceLink``, no CGNS). ``MultiZoneInterfaceSplitE2ETest`` seeds maps with
``RegisterLocalInterfaceFaceOnMaps`` and calls that pipeline when **both**
zones rebuild after one shared parent face splits:

#. Seed two zones / two shared global faces (``RegisterLocalInterfaceFaceOnMaps``).
#. Pipeline: ``InitNewInterfaceMappingBuffers`` copies inverse maps, clears ``l2gNew``.
#. ``RebuildZoneL2gNewAfterInterfaceSplit`` per zone (parent 0 -> 2,3; parent 1 kept).
#. ``BuildInverseInterfaceMapForZone`` for each zone (max two participants).
#. ``ExpandBcRecordForSplitInterfaces`` expands BC slots by child counts.
#. ``CommitInterfaceMappingBuffers`` swaps validated ``*New`` tables into the result.

**Algorithm note (orphan parents):** ``InitNewInterfaceMappingBuffers`` copies
the live inverse maps. Rebuild only touches globals present in ``l2gNew``.
A fully replaced parent global id therefore remains as a row that no ``l2g``
entry references. That is legacy behavior; changing it needs a dedicated
migration and tests, not a drive-by cleanup.

**Why compose tests?** Unit tests answer "does this function obey its contract?".
A thin multi-zone composition answers "do the contracts still line up in
production order?". That is the cheapest e2e that still avoids CGNS.


Workbench host and OneFLOW GUI
--------------------------------

The Geometry Contracts dialog is a contract-exploration surface, not a test
launcher. Both input and expected tables can be edited as JSON and restored to
the scenario defaults. The current backend is fixture-backed: it parses the
edited input but returns a known mock table result rather than recomputing the
geometry algorithm from that input. The expected tables are compared against
that fixture, and the topology sketch and optional PNG export are derived from
the mock output. This lets the scenario-editing and comparison loop develop
before binding the native geometry pipeline. A PASS here validates the edited
expected contract against the fixture only; it does not validate production C++
geometry behavior. A regression test explicitly locks the fixture-backed
boundary so the UI cannot accidentally imply that input edits trigger a real
geometry calculation.

``tools/geometry_workbench`` holds stdlib scenario JSON and a stub pipeline.
The OneFLOW Workbench home screen exposes **Geometry Contracts...**
(``GeometryContractDialog``), which loads those scenarios and runs the stub.
This is intentionally separate from **Tests...** (GoogleTest / regression).

When ``RunSymmetricInterfaceSplitMappingPipeline`` is bound, only the workbench
backend should change; the GUI dialog keeps the same table workflow.

Graphics are **optional**: CLI ``--plot-file`` and GUI **Export PNG...** write a
static figure (matplotlib Agg). Default path is text-only; no interactive plot
window is required.



Line interface split (algorithm stages)
---------------------------------------

``FaceSearch::SplitLine`` / ``TrySplitLineAtMidpoint`` is the production path
that still touches geometry services. The **table** and **arithmetic** steps are
already pure and testable without CGNS (``LineSplitPureStagesTest``).

**Flow (pure = open box, hard edge = filled box):**

.. code-block:: text

   [1 sole-zone] --> [2 midpoint] --> [3 locator FindPoint]
                                           | miss -> abort
                                           v
                              [4 child edge candidates]
                                           |
                                           v
                              [5 resolve registered faces]
                                           | empty -> abort
                                           v
                              [6 zone nodeMesh find/add]
                                           |
                                           v
                              [7 local line mapping]

**Midpoint formula** (``ComputeSegmentMidpoint``; component-wise):

.. math::

   \mathbf{m} = \frac{1}{2}(\mathbf{a} + \mathbf{b})

No tolerance is applied at this stage. Matching the global mid-node is a
separate locator query (stage 3), not part of the arithmetic.

**Stage list:**

#. ``RequireSoleZoneIdForGlobalFace`` — exactly one zone on the global face.
#. ``ComputeSegmentMidpoint`` — formula above; no mesh.
#. ``PointLocator::FindPoint`` — global mid-node must **already** exist; else abort.
#. ``BuildGlobalLineSplitCandidates`` — directed edges :math:`(p_0,m)` and :math:`(m,p_1)`.
#. ``ResolveRegisteredChildFaces`` — empty resolution aborts before mesh mutation.
#. Zone ``nodeMesh`` find-or-add (local mid index).
#. ``BuildLocalLineSplitMapping`` — parent-local rings and endpoint flags.

Header comments on the pure stages use ``@doc-algorithm line-interface-split``
and ``@doc-stage N`` so a future extractor can rebuild this section from
``FaceSearch.h`` without inventing a second algorithm description.

The remaining hard edge is intentional: stages 3 and 6 need a live locator and a
zone mesh. Do not fold mesh mutation into the pure candidate builders.



Documentation from code comments
--------------------------------

Prefer **one source of truth** for algorithm stages:

* C++ public headers carry short, structured tags (ASCII only), for example::

    // @doc-algorithm line-interface-split
    // @doc-stage 2
    // @doc-pure
    // @doc-formula m = (1/2) * (a + b) component-wise in R^3

* ``doc/source/*.rst`` holds the narrative: motivation, flow diagrams, math,
  open pressures, and ROI priority. Diagrams and equations belong here (or in
  generated pages), not in long ``//`` essays inside ``.cpp`` files.

* Tests name the contract (``LineSplitPureStagesTest``). Do not duplicate the
  full stage list in three places; link the test from the rst page.

A full auto-generator is optional. Until it exists, keep tags stable so a script
or an AI pass can collect ``@doc-algorithm`` groups into a page. Do not enable a
tool that rewrites production code comments from prose automatically.


Open design pressures
---------------------

* Quad/child-node metadata is indexed by global child face ID, but a flag of
  ``1`` means the stored node index is relative to a parent face ring. The
  consumer path is ``FaceSearch::SplitQuad2Tri`` ->
  ``FaceTopo::SetNewFace2Node`` ->
  ``RebuildFaceNodeListsAfterInterfaceSplit`` ->
  ``ResolveChildFaceNodeIds``. The final resolver receives the current parent
  ring, while its mapping row is selected by global child ID. A second parent
  that shares that child can overwrite the row with different parent-local
  indices.

  Parent-specific tables address this; the fixture is
  ``TwoParentSharedChildMappingTest`` (rings ``[10,20,30,40]`` and
  ``[20,30,40,50]``, shared child resolved to ``[20,30,40]`` under both
  parents). Merely rotating one parent's ring does not create a second global
  parent face in the registry.

  Each parent child list must contain a global child face ID at most once.
  The mapping rows are indexed by child position, so duplicate IDs would make a
  single ``SetChildMapping(parent, child, ...)`` call ambiguous: the same child ID
  would name multiple positions. ``SetChildMapping`` rejects duplicate
  registrations before resizing or writing either parent-specific row. The
  ``SplitQuad2Tri`` and ``SplitLine`` share the pure resolution step
  ``ResolveRegisteredChildFaces``: it maps candidate node rings to registered
  global face IDs, skips candidates absent from the registry, and rejects a
  repeated child ID before returning either output. The caller then publishes
  the resolved child IDs and their parent-relative mapping rows. This keeps
  candidate resolution separate from parent mapping publication and prevents
  the quad and line paths from drifting in duplicate handling. A malformed
  parent therefore has an empty rebuilt mapping after a duplicate-candidate
  failure, rather than a prefix that looks complete
  (**test:** ``FaceSearchBasicTest.QuadSplitRejectsDuplicateChildBeforeAppendingIt``).
  ``SetChildMapping`` remains a second, API-level guard for callers that register
  mappings directly (**test:** ``FaceSearchBasicTest.SetChildMappingRejectsDuplicateChildRegistration``).

  Re-registering the same parent-child mapping is idempotent, but writing a
  different local-corner row for that same relation is rejected before either
  stored row changes (**test:**
  ``FaceSearchBasicTest.SetChildMappingRejectsConflictingMappingForSameParent``).
  This is distinct from two different parents sharing one global child: those
  rows are expected to differ because each row is interpreted against its own
  parent ring. The invariant is therefore per parent-child relation, not per
  global child ID.

  Mapping collection prefers the parent-specific row whenever it contains
  data. The legacy child-global row is only a compatibility fallback when both
  parent-specific rows for that parent are empty. This distinction matters:
  a malformed or incomplete parent-specific row is an error, not a reason to
  fall back to metadata that may have been written for another parent. The
  regression ``RebuildFaceNodeListsTest.ParentSpecificMappingTakesPrecedenceOverLegacyRows``
  deliberately gives the legacy and parent-specific rows conflicting values
  and verifies that the parent-specific mapping wins.

  The rebuild consumer follows the same publication rule. Both overloads of
  ``RebuildFaceNodeListsAfterInterfaceSplit`` assemble the replacement face
  lists in a local buffer and swap it into the caller output only after every
  child mapping resolves. If a later child has an invalid parent-local corner
  or malformed row, the caller's previous output remains intact
  (test: ``RebuildFaceNodeListsTest.FailedRebuildDoesNotPublishPartialFaceLists``).

  The overload that accepts child-global node tables is now an adapter: it
  translates each interface's child IDs into the parent-local row layout and
  delegates reconstruction to the same core overload. This keeps compatibility
  at the input boundary instead of maintaining two copies of face-list
  rebuilding logic. The regression
  ``RebuildFaceNodeListsTest.LegacyMappingAdapterUsesCanonicalRebuildAtomically``
  verifies that a failure in a later legacy-mapped child does not publish a
  partial face list.

  The split table now retains parent-specific child node IDs and flags in
  ``childFaceNodeIdsByParent`` and ``childFaceNodeFlagsByParent``. Global
  child IDs remain the identity and deduplication key; each parent owns an
  aligned row for each child. ``SetChildMapping`` permits a row to lag the child
  list while children are registered incrementally, but the ID and flag rows
  must always have equal lengths and must not exceed the parent's child count.
  The legacy global child-ID and flag rows must also be length-aligned before
  a parent-specific row can be written; malformed compatibility data is
  rejected without changing either mapping domain. This prevents a successful
  parent-specific write from masking a broken global fallback row.
  It rejects a mismatched pair before resizing either row, so malformed metadata
  is not silently normalized. This invariant is covered by
  ``FaceSearchBasicTest.SetChildMappingRejectsMisalignedParentRowsWithoutMutation``. ``FaceTopo::SetNewFace2Node`` collects these
  rows and passes them to the parent-aware overload of
  ``RebuildFaceNodeListsAfterInterfaceSplit``. The regression
  ``RebuildFaceNodeListsTest.SharedChildUsesEachParentLocalMapping`` checks
  that two parents sharing child ID 7 rebuild to the same ordered global node
  list despite different parent-local indices. The older global-indexed rows
  remain temporarily for compatibility with existing helper callers. Both the
  legacy node-ID row and its flag row retain the first populated mapping for a
  shared child; later parents may have different local-corner IDs and flags but
  must not overwrite those compatibility values. New parent-aware consumers
  must use the parent-specific tables. The regression
  ``FaceSearchBasicTest.SharedChildKeepsParentSpecificMappingRows`` checks
  both the differing parent rows and the unchanged legacy rows.

  The split-table lifecycle follows the same ownership boundary. ``ResetForFaceCount``
  sizes and clears both global child tables and parent-specific rows; a reset must
  not leave mapping rows from a previous face registry. ``ClearParentMapping``
  removes only the selected parent's child list and parent-relative rows. It must
  preserve the global child metadata and any other parent's row, because the child
  can be shared. ``FaceSearch::Clear`` is broader: it invalidates the face registry
  and therefore clears every split table. The regression tests
  ``FaceSearchBasicTest.ResetForFaceCountClearsAllSplitTables`` and
  ``FaceSearchBasicTest.ClearParentMappingOnlyClearsRequestedParent`` encode these
  distinct lifecycle scopes.

  The end-to-end regression
  ``RebuildFaceNodeListsTest.SplitProducerCollectorAndRebuildPreserveSharedChildMapping``
  drives the real ``SplitQuad2Tri`` producer through
  ``CollectBoundaryInterfaceChildFaceIds``,
  ``CollectBoundaryInterfaceChildNodeMappings``, and the parent-aware rebuild
  consumer. Two distinct quad parents share one registered triangle. The first
  parent ring is ``[100, 101, 102, 103]``; the second is
  ``[101, 102, 103, 104]``. Their shared child is ``[101, 102, 103]``, with
  parent-local rows ``[1, 2, 3]`` and ``[0, 1, 2]``. Both reconstructed
  outputs must preserve the registered child ring ``[101, 102, 103]``. This
  locks child identity, parent-local mapping, interface-local collection order,
  and final node reconstruction together rather than testing only isolated
  helpers.


These remain harder to unit-test and mark follow-up work:

#. ``IFaceLink`` dual-map fields are private. Seed and observe through
   ``EnsureInterfaceZoneCount`` / ``EnsureLocalInterfaceFaceCount`` /
   ``RegisterLocalInterfaceFace`` and the Get*/Mutable* accessors. Geometry
   call sites no longer touch ``l2g`` / ``gI2Zid`` / ``g2l`` fields.
#. Line mid-point split still needs a ready locator for production grids.
   ``point_search`` is private (``GetPointLocator`` /
   ``EnsurePointLocatorReady``). Tests inject a double via
   ``ReplacePointLocator(locator, markReady=true)`` so mid-point / id
   resolution can run without nodeMesh boxes (**test:**
   ``IFaceLinkPointLocatorInjectionTest``). Zone policy remains
   ``RequireSoleZoneIdForGlobalFace`` (**test:**
   ``RequireSoleZoneIdForGlobalFaceTest``).
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
       ``ValidateBoundaryBcLayout`` /
       ``ReorderBoundaryFirstAndValidate``
     - ``FaceTopoReorderLinkTest``,
       ``ReorderBoundaryFirstAndValidateTest``
     - `Boundary-first reorder algorithm`
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
   * - ``RegisterLocalInterfaceFaceOnMaps`` / CreateLink mapping
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

Registered child-face resolution
--------------------------------

``ResolveRegisteredChildFaces`` is a lookup stage, not a mesh-generation
stage. Its caller supplies candidate node rings; the resolver asks
``FaceSearchBasic::FindFace`` which candidates already have global face IDs.
It returns pairs of ``(globalFaceId, candidateIndex)`` for matches only.

The algorithm has four important properties:

#. **Identity is order-independent.** The face lookup keys a face by its node
   set, while ``faceNodeRings`` preserves the first registered node order.
   Reversing a candidate ring therefore does not create a new global face.
#. **Missing candidates are skipped.** Resolution does not register new faces
   and does not invent IDs for candidates absent from the registry.
#. **Candidate positions are preserved.** The returned candidate index points
   back to the original candidate list. Parent-local node IDs and flags must
   use this index; compacting candidates before applying local mappings would
   silently attach the wrong mapping after a partial match.
#. **Duplicate resolved IDs are rejected.** If two candidates resolve to the
   same global face ID, the resolver throws instead of publishing an ambiguous
   parent-to-child relation.

The resolver checks duplicate resolved IDs with a hash set while keeping the
result vector in candidate order. For `k` candidate faces, this changes the
duplicate-check portion from a repeated linear scan (worst-case quadratic in
the number of resolved children) to expected linear time, without changing
the registry lookup or output ordering. The vector remains the authoritative
ordered result; the set is only an identity guard. This is an implementation
detail, not a new topology rule.

The registered-child contract is covered by
``FaceSearchBasicTest.RegisteredChildResolutionReturnsOnlyKnownChildren``,
with companion tests for partial matches, no matches, duplicate identity, and
the original candidate index. The scenario
``tools/geometry_workbench/scenarios/registered_child_resolution.json`` makes
one of these cases inspectable as input and expected tables. Its ``result``
must remain unrecorded until a native test runner captures execution; a stub
that copies expected values is not evidence that the C++ algorithm passed.

This seam is architecturally useful because it separates three decisions:
candidate topology is constructed by the split algorithm, global face identity
is resolved by the registry, and parent-relative mapping rows are published
by ``SetChildMapping``. Keeping those decisions distinct lets tests pinpoint
whether a failure comes from candidate construction, face lookup, or mapping
publication without building a complete CGNS case.
Interface split pipeline and publication boundary
--------------------------------------------------

``RunSymmetricInterfaceSplitMappingPipeline`` is a table-level composition
of the production mapping steps. It is deliberately independent of ``IFaceLink``,
CGNS objects, and mesh geometry, so its mapping contract can be tested with
small tables:

#. Copy the live maps into candidate buffers (``InitNewInterfaceMappingBuffers``).
#. Rebuild each zone's local-to-global row by replacing split parents with their
   ordered child global IDs.
#. Rebuild the inverse global-to-zone and global-to-local tables from those rows.
#. Expand boundary-condition records using the child counts.
#. Commit the candidate maps and publish the result.

The key invariant is that these tables form one logical result even though they
are stored separately. Publishing ``nChild`` or expanded BC arrays before a later
stage succeeds can expose a mixed state: new counts alongside old maps. The
pipeline now builds an ``InterfaceSplitMappingPipelineResult`` candidate and
moves it to the caller's output only after every stage succeeds. A regression
test deliberately supplies too few child counts and verifies that all preexisting
output fields remain unchanged.

This is a strong exception guarantee for the **table-level result**; it does
not make the live production ``IFaceLink`` rebuild transactional by itself.
The production path still owns its staging/commit lifecycle, and that lifecycle
must be reviewed separately when it changes.

The multi-zone end-to-end table test also exposes an important distinction:
after replacing a parent face with children, the copied inverse tables may retain
an unreferenced row for the retired parent global ID. Detection is implemented by
``CollectUnreferencedGlobalInterfaceFaces``. Detection is not deletion: global
IDs are indices into several related tables, so removing or compacting a row
would require consistently rewriting every local-to-global map and any
parent/child references. Until that full ID-remapping contract is defined, keep
the diagnostic and do not silently compact the maps.

Point-location tolerance and search stages
-------------------------------------------

The point locator uses the ADT as a broad-phase spatial index, then applies an
exact squared-distance test as the narrow phase. The query first searches an
axis-aligned box expanded by ``tolerance`` in each coordinate. Every candidate
is then tested against ``dx*dx + dy*dy + dz*dz < tolerance*tolerance``; the
first-stage box is therefore only a candidate filter, not the final matching
criterion. A point outside the exact radius is not a match even if it lies in
the search box.

For a single grid, the standard ADT setup derives tolerance as one quarter of
the minimum grid spacing. For multiple grids, it derives tolerance as one tenth
of the minimum spacing and applies the existing small-value fallback when that
result is below ``1.0e-10``. ``CreateStandardADTByTolerance`` instead uses the
caller-provided tolerance. These are distinct policies and should not be
unified without checking interface-node matching and midpoint-splitting call
sites: changing tolerance changes which geometric points are considered the
same point.

The bounding-box overload of ``PointLocator::Initialize`` only reads its input
bounds and copies them before expansion. Its inputs and the corresponding
``CreateStandardADT`` overload now use ``const RealField &`` to express that
contract; ``ShiftMinMaxBox`` remains mutable because it expands the copied
bounds in place.

The point locator's observation-only methods also expose const-qualified
contracts: ``GetNPoint``, ``GetTol``, ``GetPoint``, and ``GetFaceCoorList`` do
not modify locator state. This lets consumers inspect a locator through a
``const PointLocator &`` while leaving ``FindPoint`` and ``AddPoint`` mutable:
search updates the last-result ID, and insertion changes the point store and
spatial index. The distinction is behavioral, not cosmetic, and the injected
locator test exercises the query methods through a const base reference.

The coordinate-vector overloads of ``FindPoint`` and ``AddPoint`` now accept
``const RealField &`` as well. Both operations read the supplied coordinates;
their state changes belong to the locator (the result ID, point arrays, and ADT),
not to the caller's coordinate buffer. This avoids requiring a mutable vector
for an input-only role and makes the ownership boundary clearer.

Const-correct access and pointer ownership
------------------------------------------

The interface-link layer now provides both mutable and read-only grid access.
``IFaceLink::GetGrid`` returns ``Grid &`` on a mutable object and
``const Grid &`` on a const object. This makes observation-only code express
its intent in the type system without introducing a second grid lookup path.
The accessor is deliberately a reference: a zone index that passes the range
check must identify an owned grid, not an optional result.

The same distinction applies to collections of grid pointers. ``GridViews`` is
a non-owning mutable view; making the vector itself const does not make its
pointed-to grids const. ``ConstGridViews`` and ``AsConstGridViews`` provide a
separate read-only view whose elements are ``const Grid *``. This keeps the
legacy mutable path intact while allowing read-only algorithms to advertise
their actual contract in their parameter types. Neither view owns the grids;
the originating ``Grids`` collection must outlive every view.

For the broader const/raw-pointer audit, use semantic ownership rather than
a mechanical pointer-to-reference rewrite:

* Use ``const T &`` for required read-only inputs and ``T &`` for required
  mutable inputs.
* Add paired const/non-const accessors when callers may observe an owned object
  through either a mutable or const owner.
* Use ``std::unique_ptr<T>`` for exclusive ownership and transfer ownership
  explicitly; use ``std::shared_ptr<T>`` only when shared lifetime is real.
* Keep a raw ``T *`` when null is a meaningful state or when the pointer is a
  non-owning view whose lifetime is controlled elsewhere. Document that
  contract and check null before dereferencing.
* Do not remove a pointer merely because it is syntactically raw. In particular,
  legacy grid lookup APIs can encode missing or non-local zones with null, while
  ``IFaceLink::GetGrid`` is a required indexed access and is better expressed
  as a reference.

The purpose of this audit is to make lifetime, nullability, and mutation visible
at API boundaries. It is not a project-wide pointer count reduction exercise.
Prefer small, compile-checked changes that improve a real call path, with a
focused test for the type-level contract. The next audit targets are public
geometry accessors and grid/mesh lookup APIs; each candidate must be classified
as required, optional, owning, or non-owning before changing its type.
