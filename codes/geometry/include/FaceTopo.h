/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


#pragma once
#include "Constant.h"
#include "HXDefine.h"
#include <vector>
#include <string>
#include <fstream>
#include <memory>


BeginNameSpace( ONEFLOW )

class BcRecord;
class IFaceLink;
struct InterfaceFaceSplitTable;
class InterFace;
class BcManager;
class Grid;

class FaceTopo
{
public:
    FaceTopo();
    ~FaceTopo();
public:
    int nCells;
    IntField fTypes;
    LinkField faces;
    LinkField c2f;

public:
    IntField lPosition, rPosition;
private:
    IntField lCells, rCells;
    IntField faceFlags;
private:
    std::unique_ptr< BcManager > bcManager;
    Grid * grid = nullptr;

    // Staging buffers for topology rebuild. Empty after a successful commit
    // via UpdateOtherTopologyTerm(); non-empty only while a rebuild is pending.
    LinkField rebuildFaces;
    IntField rebuildLeftCells;
    IntField rebuildRightCells;

    void ClearRebuildBuffers();

public:
    // Preferred name for the total face count.
    HXSize_t GetFaceCount() const { return fTypes.size(); }
    // Compatibility wrapper; prefer GetFaceCount() for new code.
    HXSize_t GetNFaces() const { return GetFaceCount(); }
    HXSize_t CalcTotalFaceNodes() const;
    // Preferred name for the boundary face count.
    HXSize_t GetBoundaryFaceCount() const;
    // Compatibility wrapper; prefer GetBoundaryFaceCount() for new code.
    HXSize_t GetNBFaces() const;
    // Preferred name for setting the boundary face count.
    void SetBoundaryFaceCount( HXSize_t boundaryFaceCount );
    // Compatibility wrapper; prefer SetBoundaryFaceCount() for new code.
    void SetNBFaces( HXSize_t nBFaces );
    void BindGrid( Grid & grid );
    Grid & GetGrid();
    const Grid & GetGrid() const;
    BcRecord & GetBcRecord();
    const BcRecord & GetBcRecord() const;
    IntField & GetFaceFlags();
    const IntField & GetFaceFlags() const;
    IntField & GetLeftCells();
    const IntField & GetLeftCells() const;
    IntField & GetRightCells();
    const IntField & GetRightCells() const;
    LinkField & GetFaces();
    const LinkField & GetFaces() const;
    IntField & GetFaceTypes();
    const IntField & GetFaceTypes() const;
    void PrepareBoundaryConditions();
    bool HasInterfaceBoundary() const;
    // Grow left/right cells, face-local positions, and faceFlags to match faces.size().
    // Prefer StageNewFace on FaceSolver for full parallel growth (BC slots too).
    void EnsureFaceSideArraysSized();
    // Compatibility wrapper; prefer EnsureFaceSideArraysSized().
    void ResizeAll() { EnsureFaceSideArraysSized(); }
public:
    void ModifyFaceNodeId( IFaceLink & iFaceLink );
    void SetNewFace2Node( IFaceLink & iFaceLink );
    void SetNewFace2Cell( IFaceLink & iFaceLink );
    void ModifyBoundaryInformation( IFaceLink & iFaceLink );
    void ResetNumberOfBoundaryCondition( IFaceLink & iFaceLink );
    void ConstructNewInterfaceMap( IFaceLink & iFaceLink );
    void UpdateOtherTopologyTerm();
    // True while rebuildFaces / rebuildLeftCells / rebuildRightCells hold
    // staged data that has not yet been committed (or was not cleared).
    bool HasPendingTopologyRebuild() const;
    void GenerateI2B( InterFace & interFace );
public:
    bool GetSId( int iFace, int iPosition, int & sId );
    bool GetTId( int iFace, int iPosition, int & tId );
    void CalcC2C( LinkField & c2c );
    void ReorderLink();
};

// Build newToOld permutation: faces with right cell INVALID (boundary) first,
// then interior faces. Relative order within each group is preserved.
// newToOld[newIndex] == oldIndex.
IntField BuildBoundaryFirstPermutation( const IntField & rightCells );

// Count faces whose right cell is INVALID_INDEX.
HXSize_t CountBoundaryFaces( const IntField & rightCells );

// After ReorderLink + AppendBoundaryBcToRecord: first nBFaces faces must be
// boundary, the rest interior, and BcRecord length must equal nBFaces.
// Throws std::logic_error when the layout is inconsistent.
void ValidateBoundaryBcLayout( const FaceTopo & faceTopo );

// One post-stage transaction: boundary-first ReorderLink then
// ValidateBoundaryBcLayout. Call after BcRecord has been filled from staged
// faces in discovery order of boundary faces.
void ReorderBoundaryFirstAndValidate( FaceTopo & faceTopo );


// Expand boundary BC arrays when interface faces may split into children.
// childCountPerInterfaceLocal is parallel to interface faces in boundary-face
// order (same order as walking bcType with IsInterfaceBc). For each interface
// entry, emit max(childCount, 1) copies of that BC; non-interface entries emit once.
void ExpandBcRecordForSplitInterfaces(
    const IntField & bcType,
    const IntField & bcNameId,
    const IntField & childCountPerInterfaceLocal,
    IntField & outBcType,
    IntField & outBcNameId );

// Rebuild one zone's inverse references from l2gNewForZone. Existing slots
// for zoneId are refreshed with the new local index; absent zones are appended
// only while the legacy two-participant capacity has room.
void BuildInverseInterfaceMapForZone(
    int zoneId,
    const IntField & l2gNewForZone,
    LinkField & gI2ZidNew,
    LinkField & g2lNew );

// Rebuild one zone's l2gNew and nChild after interface face splits.
// parentGlobalFaceIds[i] is the global id of parent interface local i;
// childGlobalIdsPerParent[i] is empty (keep parent) or the replacement child
// global ids. nChildForZone starts as one slot per parent (child count or 0)
// and appends a 0 for each emitted child entry (legacy layout).
void RebuildZoneL2gNewAfterInterfaceSplit(
    const IntField & parentGlobalFaceIds,
    const LinkField & childGlobalIdsPerParent,
    IntField & l2gNewForZone,
    IntField & nChildForZone );

// Rebuild left/right cell arrays after interface boundary faces may split.
// For boundary faces (0..nBFaces-1): interface faces emit max(childCount,1)
// copies of the original left cell; other boundary faces emit once. New
// boundary right cells are ghost indices iFace + nCells. Interior faces
// (nBFaces..nFaces-1) copy left/right unchanged.
void RebuildFaceCellSidesAfterInterfaceSplit(
    const IntField & bcType,
    const IntField & leftCells,
    const IntField & rightCells,
    int nBFaces,
    int nFaces,
    const IntField & childCountPerInterfaceLocal,
    int nCells,
    IntField & outLeftCells,
    IntField & outRightCells );

// Walk boundary bcType[0..nBFaces) and, for each interface entry, append the
// child global-id list from GetInterfaceChildFaceIds (interface-local order).
void CollectBoundaryInterfaceChildFaceIds(
    const IntField & bcType,
    int nBFaces,
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    LinkField & outChildGlobalIdsPerInterfaceLocal );

// Derived from CollectBoundaryInterfaceChildFaceIds (sizes only).
void CollectBoundaryInterfaceChildCounts(
    const IntField & bcType,
    int nBFaces,
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    IntField & outChildCountPerInterfaceLocal );

void CollectBoundaryInterfaceChildNodeMappings(
    const IntField & bcType,
    int nBFaces,
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    std::vector< LinkField > & outNodeIdsPerInterfaceLocal,
    std::vector< LinkField > & outNodeFlagsPerInterfaceLocal );

// Rebuild face-to-node lists after interface boundary faces may split.
// childGlobalIdsPerInterfaceLocal is parallel to interface faces in boundary
// order; empty means keep the parent node list, non-empty means emit one
// resolved list per child (via ResolveChildFaceNodeIds). Interior faces are
// copied unchanged.
void RebuildFaceNodeListsAfterInterfaceSplit(
    const LinkField & faces,
    const IntField & bcType,
    int nBFaces,
    int nFaces,
    const LinkField & childGlobalIdsPerInterfaceLocal,
    const LinkField & childFaceNodeIds,
    const LinkField & childFaceNodeFlags,
    bool hasNodeMesh,
    HXSize_t nodeMeshCount,
    LinkField & outFaces );

void RebuildFaceNodeListsAfterInterfaceSplit(
    const LinkField & faces,
    const IntField & bcType,
    int nBFaces,
    int nFaces,
    const LinkField & childGlobalIdsPerInterfaceLocal,
    const std::vector< LinkField > & childNodeIdsPerInterfaceLocal,
    const std::vector< LinkField > & childNodeFlagsPerInterfaceLocal,
    bool hasNodeMesh,
    HXSize_t nodeMeshCount,
    LinkField & outFaces );


// Look up child face IDs for a local interface face and validate their node tables.
// Preferred overload: depends only on l2g + InterfaceFaceSplitTable (no IFaceLink).
// local interface face -> global parent id -> childFaceIdsOfParent[parent].
// Throws if indices are out of range or child node id/flag lengths differ.
const IntField & GetInterfaceChildFaceIds(
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    int localInterfaceFaceId );

// Convenience wrapper for existing IFaceLink call sites.
const IntField & GetInterfaceChildFaceIds(
    const IFaceLink & interfaceLink,
    int zoneId,
    int localInterfaceFaceId );

// Resolve one child face's node indices from parent-face relative and absolute IDs.
// flag == 1: childNodeIds[i] is an index into parentFaceNodes.
// otherwise: childNodeIds[i] is an absolute mesh node ID (requires hasNodeMesh).
IntField ResolveChildFaceNodeIds(
    const IntField & parentFaceNodes,
    const IntField & childNodeIds,
    const IntField & childNodeFlags,
    bool hasNodeMesh,
    HXSize_t nodeMeshCount );

// ---------------------------------------------------------------------------
// Table-level interface split mapping pipeline (workbench / tests / future GUI).
// No IFaceLink, no CGNS, no mesh. Same steps as production dual-map rebuild:
// InitNew buffers -> per-zone RebuildZoneL2gNew -> BuildInverse -> ExpandBc -> Commit.
// childGlobalIdsPerParent is applied identically to every zone (symmetric case).
// ---------------------------------------------------------------------------
struct InterfaceSplitMappingPipelineResult
{
    LinkField l2g;
    LinkField gI2Zid;
    LinkField g2l;
    LinkField nChild;
    IntField outBcType;
    IntField outBcNameId;
};

// Global face row indices in [0, globalFaceRowCount) that are never referenced
// by any local-to-global entry. After a split rebuild that keeps InitNew's inverse
// copy, retired parent globals often appear here (orphan inverse rows).
// On validation failure, outUnreferencedSorted is left unchanged.
void CollectUnreferencedGlobalInterfaceFaces(
    const LinkField & localToGlobal,
    std::size_t globalFaceRowCount,
    IntField & outUnreferencedSorted );

void RunSymmetricInterfaceSplitMappingPipeline(
    const LinkField & l2gIn,
    const LinkField & gI2ZidIn,
    const LinkField & g2lIn,
    const IntField & parentGlobalFaceIds,
    const LinkField & childGlobalIdsPerParent,
    const IntField & bcType,
    const IntField & bcNameId,
    const IntField & childCountPerInterfaceLocal,
    InterfaceSplitMappingPipelineResult & out );


EndNameSpace
