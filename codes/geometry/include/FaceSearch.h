/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )

class FaceSort
{
public:
    FaceSort();
    FaceSort( const IntField & nodeId, int fId = 0 );
    ~FaceSort();
public:
    int fId;
    IntField nodeId;
    IntField sortedNodeId;
};

class CompareFace
{
public:
    bool operator() ( const FaceSort * lhs, const FaceSort * rhs ) const
    {
        return lhs->sortedNodeId < rhs->sortedNodeId;
    }
};

class FaceSearchBasic
{
public:
    FaceSearchBasic();
    ~FaceSearchBasic();
public:
    HXVector< FaceSort * > faceArray;
    set< FaceSort * , CompareFace > faceSet;
public:
    int AddFace( const IntField & faceNode );
    int FindFace( const IntField & faceNode );
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
    void SplitQuad2Tri( FaceSort * pFaceSort );
    void SplitLine( FaceSort * pFaceSort );
    void GetLocalTri( LinkField & localTriId, LinkField & localTriFlag );
    void GetTriId( FaceSort * pFaceSort, LinkField & localTriId, LinkField & triId );
    bool GetLine( FaceSort * pFaceSort, LinkField & localLineId, LinkField & localLineFlag, LinkField & lineId );
};

EndNameSpace