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

    [[nodiscard]] bool StagedBcMatchesFaceCount() const;
public:
    void ResizeAll();
    void ScanPolygonFace( CgnsSection & cgnsSection );
    void ScanPolyhedronElement( CgnsSection & cgnsSection );
private:
    std::unique_ptr< FaceTopo > faceTopo;
};

EndNameSpace
