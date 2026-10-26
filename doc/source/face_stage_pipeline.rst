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

BC classification and CGNS mapping
----------------------------------

* After export, CGNS type ids are converted with ``BcTypeMap::MapFromCgnsInPlace``
  (**test:** ``BcTypeMapTest``).
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
  leaves the face explicitly unmatched instead of retaining a stale partner.
  The periodic partner row is scanned for the first participant other than the
  current zone/local-face pair. A self-only row leaves the face unmatched, and
  mismatched zone/local row lengths are rejected (**tests:**
  ``FindPeriodicInterfacePartnerTest``).
* ``CountInterfaceFacesMissingPartner`` counts local interface faces whose
  global face has fewer than two zone slots (**test:** ``CountInterfaceFacesMissingPartnerTest``).
  The runtime diagnostic is therefore labeled "Interface faces missing a partner", not "Periodic boundary
  faces": this counter examines the interface mapping tables and does not
  independently classify periodic boundaries. Faces with more than two zone
  slots remain excluded for legacy compatibility; that policy is not changed
  by the diagnostic wording correction.
* ``ExpandBcRecordForSplitInterfaces`` repeats BC slots when an interface face
  splits into children (**test:** ``ExpandBcRecordSplitTest``).
* ``RegisterLocalInterfaceFace`` constructs the initial global-face participant
  rows. A global face supports at most two zone participants; attempting to
  register a third participant now throws before changing ``l2g``, ``gI2Zid``, or ``g2l``
  (**test:** ``RegisterLocalInterfaceFaceTest``).
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
  per-interface child global-id lists; ``CollectBoundaryInterfaceChildCounts`` is the size-only view
  (**test:** ``CollectBoundaryInterfaceChildCountsTest``).
* ``RebuildFaceNodeListsAfterInterfaceSplit`` rebuilds face-to-node lists after
  interface splits (**test:** ``RebuildFaceNodeListsTest``).
* Quad diagonal candidates:
  ``BuildQuadDiagonalLocalTris`` creates the four local triangle candidates
  for both diagonals; ``MapLocalTrisToGlobalNodes`` maps their three local
  corners through the parent face's node ring. The mapping helper validates
  every candidate row has exactly three corners and every corner index is in
  range before replacing the output, so malformed later rows cannot leave a
  partial mapping (**test:** ``QuadDiagonalSplitTest``).
* ``SplitQuad2Tri`` clears the selected parent's previous child list before
  checking its shape; a direct call on a non-quad face therefore cannot preserve
  a stale split result (**test: ``FaceSearchBasicTest.NonQuadSplitClearsStaleParentMapping``).


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
   migrated and compatibility is verified.
