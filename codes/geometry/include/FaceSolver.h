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
#include "HXDefine.h"
#include "HXCgns.h"
#include "HXLookup.h"
#include <vector>
#include <set>
#include <map>
#include <memory>


BeginNameSpace( ONEFLOW )
class MeshPointManager;
class ElemFeature;
class FaceTopo;
class CgnsSection;

class FaceSolver
{
public:
    FaceSolver();
    ~FaceSolver();
public:
    HXLookup<int> faceLookup;
    // [Refactored] Changed from raw pointers to value types.
    // They are now automatically managed by the compiler (stack allocation).
    // Staged BC attributes per face (aligned with FaceTopo faces).
    // Copied into BcRecord as bcNameId / bcType in GridElem::GenerateCalcElement.
    IntField faceBcNameIds;
    IntField faceBcTypes;
    LinkField childFid;
public:
    FaceTopo & GetFaceTopo();
    const FaceTopo & GetFaceTopo() const;
    [[nodiscard]] std::unique_ptr< FaceTopo > ReleaseFaceTopo() noexcept;
public:
    bool CheckBcFace( IntSet & bcVertex, IntField & nodeId );
    void ScanElementFace( CgIntField & eNodeId, int eType, int eId );
    void ScanBcFace( IntSet & bcVertex, int bcType, int bcNameId );
    void ScanBcFaceDetail( IntSet & bcVertex, int bcType, int bcNameId );
    void ScanInterfaceBc();
    int GetNSimpleFace();

    // Stage a newly discovered unique face: topology, face flag, and default
    // BC slots together. Keeps faceBcTypes / faceBcNameIds / faceFlags length
    // equal to FaceTopo face count. leftCellId may be INVALID_INDEX when the
    // owning cell is assigned later (polygon path).
    void StageNewFace( const IntField & faceNodes, int faceType, int leftCellId );

    // Assign BC on an already staged face (bounds-checked). Sole write path
    // for faceBcTypes / faceBcNameIds after staging.
    void SetStagedFaceBc( int faceId, int bcType, int bcNameId );

    // First cell to claim an unclaimed staged face (polygon path). Sets the
    // claimed flag, left cell, clears right cell, and resets BC defaults via
    // SetStagedFaceBc. Throws if faceId is out of range or already claimed.
    void ClaimStagedFaceAsLeft( int faceId, int cellId );

    // Second cell sharing a claimed face. Sets right cell only.
    // Throws if faceId is out of range or the face is still unclaimed.
    void AttachStagedFaceAsRight( int faceId, int cellId );

    // True when BC slots, face flags, left/right cells, and face-local
    // positions all match the FaceTopo face count.
    [[nodiscard]] bool StagedBcMatchesFaceCount() const;

    // Face flag values used by the polygon ownership path.
    static constexpr int kFaceUnclaimed = 0;
    static constexpr int kFaceClaimed = 1;
public:
    // Safety net: resize FaceTopo cell/position arrays and FaceSolver BC/child
    // slots to the current face count. Prefer StageNewFace, which grows all
    // parallel arrays together; this remains for paths that only grew faces.
    void EnsureParallelFaceArraysSized();

    // Stage unique polygon faces from CSR-style offsets/connectivity without
    // requiring a CgnsSection. Invalid spans throw before any mutation.
    // Each new unique face is staged with INVALID left cell (owner assigned
    // later via ClaimStagedFaceAsLeft).
    void StagePolygonFacesFromConnectivity(
        int nElement,
        const CgIntField & offsets,
        const CgIntField & connectivity,
        int faceType );

    // Assign polyhedron cells to already staged faces from CSR face references.
    // Validates spans, face id range, no duplicate face per cell, and at most
    // two cells per face before any mutation. Uses ClaimStagedFaceAsLeft /
    // AttachStagedFaceAsRight only.
    void ClaimPolyhedronCellsFromConnectivity(
        int nElement,
        const CgIntField & offsets,
        const CgIntField & connectivity );

    void ScanPolygonFace( CgnsSection & cgnsSection );
    void ScanPolyhedronElement( CgnsSection & cgnsSection );
private:
    std::unique_ptr< FaceTopo > faceTopo;
};

EndNameSpace
