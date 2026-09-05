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

//class FaceSort
//{
//public:
//    FaceSort();
//    FaceSort( const IntField & nodeId, int fId = 0 );
//    ~FaceSort();
//public:
//    int fId;
//    IntField nodeId;
//    IntField sortedNodeId;
//};

class FaceSearchBasic
{
public:
    FaceSearchBasic();
    ~FaceSearchBasic();
public:
    // Store original node order of each face (index == faceId)
    HXVector< IntField > faceArray;
    //HXVector< FaceSort * > faceArray;
public:
    // Add a face. Returns the unique face ID (always >= 0).
    int AddFace( const IntField & faceNode );

    // Only find. Returns INVALID_INDEX if not found.
    int FindFace( const IntField & faceNode ) const;

    std::size_t Size() const { return lookup_.Size(); }
    void Clear();

private:
    HXLookup<int> lookup_;   // responsible for uniqueness and ID assignment
};

class IFaceLink;

class FaceSearch : public FaceSearchBasic
{
public:
    FaceSearch();
    ~FaceSearch();
public:
    IntField status;
    LinkField cFaceId;
    LinkField rCNodeId;
    LinkField rCNodeFlag;
    IFaceLink * iFaceLink;
    int gFid;
public:
    void CalcNewFaceId( IFaceLink * iFaceLink );
    //void SplitQuad2Tri( FaceSort * pFaceSort );
    void SplitQuad2Tri( int faceId );
    //void SplitLine( FaceSort * pFaceSort );
    void SplitLine( int faceId );
    void GetLocalTri( LinkField & localTriId, LinkField & localTriFlag );
    //void GetTriId( FaceSort * pFaceSort, LinkField & localTriId, LinkField & triId );
    void GetTriId( const IntField & nodeId, LinkField & localTriId, LinkField & triId );
    //bool GetLine( FaceSort * pFaceSort, LinkField & localLineId, LinkField & localLineFlag, LinkField & lineId );
    bool GetLine( const IntField & nodeId, LinkField & localLineId, LinkField & localLineFlag, LinkField & lineId );
};

EndNameSpace
