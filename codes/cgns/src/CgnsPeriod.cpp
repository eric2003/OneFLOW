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

#include "CgnsPeriod.h"
#include "FaceSearch.h"
#include "HXMath.h"
#include "NodeMesh.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

F2FMap f2fmap;

F2FMap::F2FMap()
{
    pointBasic = new PointBasic();
    faceSearchBasic = new FaceSearchBasic();
}

F2FMap::~F2FMap()
{
    delete pointBasic;
    delete faceSearchBasic;
}

void F2FMap::AddFacePair(int faceId1, int faceId2)
{
    this->face_pair.insert( pair<int, int>(faceId1, faceId2) );
}

int F2FMap::FindPeriodFace( int faceId )
{
    map< int, int >::iterator iter = this->face_pair.find( faceId );
    int face_id = -1;
    if ( iter != this->face_pair.end() )
    {
        face_id = iter->second;

    }
    return face_id;
}

void F2FMap::AddFacePoint( CgIntField & fNodeId1, CgIntField & fNodeId2, NodeMesh *nodeMesh1, NodeMesh *nodeMesh2 )
{
    IntField face1, face2;
    int fId1 = this->AddFacePoint( fNodeId1, nodeMesh1, face1 );
    int fId2 = this->AddFacePoint( fNodeId2, nodeMesh2, face2 );

    this->AddFacePair( fId1, fId2 );

    this->faceList1.push_back( face1 );
    this->faceList2.push_back( face2 );
}

int F2FMap::AddFacePoint( CgIntField & fNodeId, NodeMesh *nodeMesh, IntField & newFaceNodeId )
{
    //IntField face;
    for ( int i = 0; i < fNodeId.size(); ++ i )
    {
        int ip = fNodeId[ i ];

        Real xm = nodeMesh->xN[ ip ];
        Real ym = nodeMesh->yN[ ip ];
        Real zm = nodeMesh->zN[ ip ];

        int pid = this->pointBasic->AddPoint( xm, ym, zm );
        newFaceNodeId.push_back( pid );
    }

    int fId = this->faceSearchBasic->AddFace( newFaceNodeId );
    return fId;
}

void F2FMap::FindFace( RealField &xList, RealField &yList, RealField &zList, RealField &xxList, RealField &yyList, RealField &zzList )
{
    int nNode = xList.size();
    IntField face;
    for ( int i = 0; i < nNode; ++ i )
    {
        Real xm = xList[ i ];
        Real ym = yList[ i ];
        Real zm = zList[ i ];

        int pid = this->pointBasic->FindPoint( xm, ym, zm );
        face.push_back( pid );
    }

    int face_id = this->faceSearchBasic->FindFace( face );
    int face_period = this->FindPeriodFace( face_id );
    //if ( face_period == -1 )
    FaceSort * faceSort = this->faceSearchBasic->faceArray[ face_period ];
    IntField & node_period = faceSort->nodeId;

    for ( int i = 0; i < node_period.size(); ++ i )
    {
        int  ip = node_period[ i ];
        PointBasic::PointType & pt = this->pointBasic->pointList[ ip ];

        xxList.push_back( pt.x );
        yyList.push_back( pt.y );
        zzList.push_back( pt.z );
    }
}

#endif
EndNameSpace