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
#include "PointFactory.h"
#include "HXCgns.h"
#include <map>
using namespace std;

BeginNameSpace(ONEFLOW)

#ifdef ENABLE_CGNS

class PointBasic;
class FaceSearchBasic;
class NodeMesh;

class F2FMap
{
public:
    F2FMap();
    ~F2FMap();
public:
    PointBasic * pointBasic;
    FaceSearchBasic * faceSearchBasic;
    LinkField faceList1, faceList2;
    map< int, int > face_pair;
public:
    void AddFacePair(int faceId1, int faceId2);
    void AddFacePoint( CgIntField & fNodeId1, CgIntField & fNodeId2, NodeMesh *nodeMesh1, NodeMesh *nodeMesh2 );
    int AddFacePoint( CgIntField & fNodeId, NodeMesh *nodeMesh, IntField & newFaceNodeId );
public:
    void FindFace( RealField &xList, RealField &yList, RealField &zList, RealField &xxList, RealField &yyList, RealField &zzList );
    int FindPeriodFace( int faceId );
};

extern F2FMap f2fmap;

#endif
EndNameSpace