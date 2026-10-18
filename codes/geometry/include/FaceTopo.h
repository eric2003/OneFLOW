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

    HXSize_t nBFaces;
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
    void ResizeAll();
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

// Look up child face IDs for a local interface face and validate their node tables.
// Returns a reference into interfaceLink.face_search->cFaceId[globalFaceId].
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

EndNameSpace
