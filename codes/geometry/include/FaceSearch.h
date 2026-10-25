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
#include "HXLookup.h"
#include <set>


BeginNameSpace( ONEFLOW )

// Unique face registry: identity via HXLookup (sorted nodes); first-seen
// node ring kept privately for orientation-sensitive consumers.
class FaceSearchBasic
{
public:
    FaceSearchBasic() = default;
    ~FaceSearchBasic() = default;

    // Add a face. Returns the unique face ID (always >= 0).
    // Same node set in any order yields the same id.
    int AddFace( const IntField & faceNode );

    // Only find. Returns INVALID_INDEX if not found (order-invariant).
    int FindFace( const IntField & faceNode ) const;

    // First-seen node ring for faceId (bounds-checked). Not overwritten on
    // later AddFace of the same identity.
    const IntField & GetFaceNodeRing( int faceId ) const;

    std::size_t Size() const { return lookup_.Size(); }
    void Clear();

private:
    HXLookup<int> lookup_;   // responsible for uniqueness and ID assignment
    // Original node order of each face (index == faceId).
    HXVector< IntField > faceArray;
};

class IFaceLink;

// Tables produced when interface faces are split (e.g. quad->tri).
// Separate from FaceSearchBasic (unique face registry).
struct InterfaceFaceSplitTable
{
    LinkField childFaceIdsOfParent;
    LinkField childFaceNodeIds;
    LinkField childFaceNodeFlags;

    void Clear()
    {
        childFaceIdsOfParent.clear();
        childFaceNodeIds.clear();
        childFaceNodeFlags.clear();
    }

    void ResetForParentFaceCount( std::size_t nParentFaces )
    {
        childFaceIdsOfParent.assign( nParentFaces, IntField{} );
        childFaceNodeIds.assign( nParentFaces, IntField{} );
        childFaceNodeFlags.assign( nParentFaces, IntField{} );
    }
};

class FaceSearch : public FaceSearchBasic
{
public:
    FaceSearch() = default;
    ~FaceSearch() = default;
public:
    IntField status;
    InterfaceFaceSplitTable interfaceSplit;
    IFaceLink * iFaceLink = nullptr;
    int gFid = 0;
public:
    void CalcNewFaceId( IFaceLink & iFaceLink );
    void SplitQuad2Tri( int faceId );
    void SplitLine( int faceId );
    void GetLocalTri( LinkField & localTriId, LinkField & localTriFlag );
    void GetTriId( const IntField & nodeId, LinkField & localTriId, LinkField & triId );
    bool GetLine( const IntField & nodeId, LinkField & localLineId, LinkField & localLineFlag, LinkField & lineId );
};

// Pure helpers for quad diagonal split (testable without IFaceLink).
// BuildQuadDiagonalLocalTris clears outputs then fills four candidate tris:
// diagonal 0-2 (tris 0-1-2 and 2-3-0) and diagonal 1-3 (tris 3-0-1 and 1-2-3).
void BuildQuadDiagonalLocalTris( LinkField & localTriId, LinkField & localTriFlag );

// Map local corner indices in localTriId through nodeId into global node ids.
void MapLocalTrisToGlobalNodes(
    const IntField & nodeId,
    const LinkField & localTriId,
    LinkField & triId );

EndNameSpace
