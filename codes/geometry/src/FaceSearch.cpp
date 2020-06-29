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

#include "FaceSearch.h"
#include "PointSearch.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "IFaceLink.h"
#include "Dimension.h"
#include "Stop.h"
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

FaceSort::FaceSort()
{
    ;
}
FaceSort::FaceSort( const IntField & nodeId, int fId )
{
    this->fId         = fId;
    this->nodeId       = nodeId;
    this->sortedNodeId = nodeId;
    std::sort( sortedNodeId.begin(), sortedNodeId.end() );
}

FaceSort::~FaceSort()
{
    ;
}


FaceSearchBasic::FaceSearchBasic()
{
    ;
}

FaceSearchBasic::~FaceSearchBasic()
{
    ;
}

int FaceSearchBasic::AddFace( const IntField & faceNode )
{
    int fId = faceArray.size();
    FaceSort * faceSort = new FaceSort( faceNode, fId );

    set< FaceSort *, CompareFace >::iterator iter = faceSet.find( faceSort );

    if ( iter == faceSet.end() )
    {
        faceSet.insert( faceSort );
        faceArray.push_back( faceSort );
        return fId;
    }
    else
    {
        delete faceSort;
        return ( * iter )->fId;
    }
}

int FaceSearchBasic::FindFace( const IntField & faceNode )
{
    int fId = faceArray.size();
    FaceSort * faceSort = new FaceSort( faceNode, fId );

    set< FaceSort *, CompareFace >::iterator iter = faceSet.find( faceSort );
    int face_id = -1;
    if ( iter != faceSet.end() )
    {
        face_id = ( * iter )->fId;

    }
    delete faceSort;
    return face_id;
}

FaceSearch::FaceSearch()
{
    ;
}

FaceSearch::~FaceSearch()
{
    ;
}

void FaceSearch::CalcNewFaceId( IFaceLink * iFaceLink )
{
    this->iFaceLink = iFaceLink;
    int nFace = this->faceArray.size();
    this->status.resize( nFace, -1 );
    this->cFaceId.resize( nFace );
    this->rCNodeId.resize( nFace );
    this->rCNodeFlag.resize( nFace );

    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        for ( int iFace = 0; iFace < nFace; ++ iFace )
        {
            this->gFid = iFace;
            FaceSort * faceSort = this->faceArray[ iFace ];
            this->SplitQuad2Tri( faceSort );
        }
    }
    else if ( Dim::dimension == ONEFLOW::TWO_D )
    {
        for ( int iFace = 0; iFace < nFace; ++ iFace )
        {
            this->gFid = iFace;
            FaceSort * faceSort = this->faceArray[ iFace ];
            this->SplitLine( faceSort );
        }
    }
}

void FaceSearch::SplitQuad2Tri( FaceSort * pFaceSort )
{
    int nNode = pFaceSort->nodeId.size();
    if ( nNode <= 3 ) return;

    LinkField localTriId, localTriFlag;
    this->GetLocalTri( localTriId, localTriFlag );

    LinkField triId;
    this->GetTriId( pFaceSort, localTriId, triId );

    for ( int iTri = 0; iTri < triId.size(); ++ iTri )
    {
        IntField & tri = triId[ iTri ];

        FaceSort * faceSort = new FaceSort( tri );
        set< FaceSort *, CompareFace >::iterator iter = this->faceSet.find( faceSort );

        if ( iter != this->faceSet.end() )
        {
            int fId = ( * iter )->fId;
            int pFid = pFaceSort->fId;

            this->cFaceId[ pFid ].push_back( fId );

            this->rCNodeId  [ fId ] = localTriId  [ iTri ];
            this->rCNodeFlag[ fId ] = localTriFlag[ iTri ];
        }

        delete faceSort;
    }
}

void FaceSearch::SplitLine( FaceSort * pFaceSort )
{
    int nNode = pFaceSort->nodeId.size();
    if ( nNode >= 3 ) return;

    LinkField localLineId, localLineFlag;
    LinkField lineId;

    if ( ! this->GetLine( pFaceSort, localLineId, localLineFlag, lineId ) ) return;

    for ( int iLi = 0; iLi < lineId.size(); ++ iLi )
    {
        IntField & lId = lineId[ iLi ];

        FaceSort * faceSort = new FaceSort( lId );
        set< FaceSort *, CompareFace >::iterator iter = this->faceSet.find( faceSort );

        if ( iter != this->faceSet.end() )
        {
            int fId = ( * iter )->fId;
            int pFid = pFaceSort->fId;

            this->cFaceId[ pFid ].push_back( fId );

            this->rCNodeId  [ fId ] = localLineId  [ iLi ];
            this->rCNodeFlag[ fId ] = localLineFlag[ iLi ];
        }

        delete faceSort;
    }
}

void FaceSearch::GetLocalTri( LinkField & localTriId, LinkField & localTriFlag )
{
    IntField localId1( 3 ), localId2( 3 );
    localId1[ 0 ] = 0;
    localId1[ 1 ] = 1;
    localId1[ 2 ] = 2;

    localId2[ 0 ] = 2;
    localId2[ 1 ] = 3;
    localId2[ 2 ] = 0;

    IntField localFlag1( 3 ), localFlag2( 3 );

    localFlag1[ 0 ] = 1;
    localFlag1[ 1 ] = 1;
    localFlag1[ 2 ] = 1;

    localFlag2[ 0 ] = 1;
    localFlag2[ 1 ] = 1;
    localFlag2[ 2 ] = 1;

    localTriId.push_back( localId1 );
    localTriId.push_back( localId2 );

    localTriFlag.push_back( localFlag1 );
    localTriFlag.push_back( localFlag2 );

    //The second splitting method
    localId1[ 0 ] = 3;
    localId1[ 1 ] = 0;
    localId1[ 2 ] = 1;

    localId2[ 0 ] = 1;
    localId2[ 1 ] = 2;
    localId2[ 2 ] = 3;

    localFlag1[ 0 ] = 1;
    localFlag1[ 1 ] = 1;
    localFlag1[ 2 ] = 1;

    localFlag2[ 0 ] = 1;
    localFlag2[ 1 ] = 1;
    localFlag2[ 2 ] = 1;

    localTriId.push_back( localId1 );
    localTriId.push_back( localId2 );

    localTriFlag.push_back( localFlag1 );
    localTriFlag.push_back( localFlag2 );

}

void FaceSearch::GetTriId( FaceSort * pFaceSort, LinkField & localTriId, LinkField & triId )
{
    triId.resize( localTriId.size() );

    for ( int iTri = 0; iTri < localTriId.size(); ++ iTri )
    {
        for ( int iNode = 0; iNode < 3; ++ iNode )
        {
            triId[ iTri ].push_back( pFaceSort->nodeId[ localTriId[ iTri ][ iNode ] ] );
        }
    }
}

bool FaceSearch::GetLine( FaceSort * pFaceSort, LinkField & localLineId, LinkField & localLineFlag, LinkField & lineId )
{
    PointSearch * point_search = this->iFaceLink->point_search;

    int p0 = pFaceSort->nodeId[ 0 ];
    int p1 = pFaceSort->nodeId[ 1 ];

    RealField coor0( 3 ), coor1( 3 ), coor2( 3 );

    point_search->GetPoint( p0, coor0[ 0 ], coor0[ 1 ], coor0[ 2 ] );
    point_search->GetPoint( p1, coor1[ 0 ], coor1[ 1 ], coor1[ 2 ] );

    for ( int m = 0; m < 3; ++ m )
    {
        coor2[ m ] = half * ( coor0[ m ] + coor1[ m ] );
    }

    if ( point_search->FindPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] ) == -1 ) return false;

    int nNZone = this->iFaceLink->gI2Zid[ this->gFid ].size();

    if ( nNZone > 1 )
    {
        Stop( "impossible" );
    }

    int zoneIndex = this->iFaceLink->gI2Zid[ this->gFid ][ 0 ];
    Grid * grid = ( this->iFaceLink->grids )[ zoneIndex ];

    int nNode = grid->nodeMesh->GetNumberOfNodes();
    int pId = nNode;

    grid->nodeMesh->AddPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] );

    int p2 = point_search->AddPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] );

    IntField id1( 2 ), id2( 2 );
    IntField localId1( 2 ), localId2( 2 );
    IntField localFlag1( 2 ), localFlag2( 2 );

    localId1[ 0 ] = 0;
    localId1[ 1 ] = pId;

    localId2[ 0 ] = pId;
    localId2[ 1 ] = 1;

    localLineId.push_back( localId1 );
    localLineId.push_back( localId2 );

    localFlag1[ 0 ] =   1;
    localFlag1[ 1 ] = - 1;

    localFlag2[ 0 ] = - 1;
    localFlag2[ 1 ] =   1;

    localLineFlag.push_back( localFlag1 );
    localLineFlag.push_back( localFlag2 );

    id1[ 0 ] = p0;
    id1[ 1 ] = p2;

    id2[ 0 ] = p2;
    id2[ 1 ] = p1;

    lineId.push_back( id1 );
    lineId.push_back( id2 );
    return true;
}

EndNameSpace