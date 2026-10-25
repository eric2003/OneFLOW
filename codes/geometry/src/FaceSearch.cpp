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

#include "FaceSearch.h"
#include "PointLocator.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "IFaceLink.h"
#include "Dimension.h"
#include "Fatal.h"
#include <algorithm>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

int FaceSearchBasic::AddFace( const IntField & faceNode )
{
    // HXLookup automatically sorts the nodes for uniqueness
    auto [faceId, isNew] = lookup_.FindOrAdd( faceNode );

    if ( isNew )
    {
        // New face: store the original node order
        faceArray.push_back( faceNode );
    }

    return faceId;   // always return a valid ID
}

int FaceSearchBasic::FindFace( const IntField & faceNode ) const
{
    return lookup_.Find( faceNode );   // returns INVALID_INDEX if not found
}

void FaceSearchBasic::Clear()
{
    faceArray.clear();
    lookup_.Clear();
}

const IntField & FaceSearchBasic::GetFaceNodeRing( int faceId ) const
{
    if ( faceId < 0 || static_cast< size_t >( faceId ) >= this->faceArray.size() )
    {
        throw std::out_of_range( "FaceSearchBasic::GetFaceNodeRing: face ID is out of range" );
    }
    return this->faceArray[ static_cast< size_t >( faceId ) ];
}


void FaceSearch::CalcNewFaceId( IFaceLink & iFaceLink )
{
    this->iFaceLink = &iFaceLink;
    int nFaces = static_cast< int >( this->Size() );
    this->status.assign( nFaces, -1 );
    this->interfaceSplit.ResetForParentFaceCount( static_cast< std::size_t >( nFaces ) );

    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        for ( int iFace = 0; iFace < nFaces; ++ iFace )
        {
            this->gFid = iFace;
            //FaceSort * faceSort = this->faceArray[ iFace ];
            this->SplitQuad2Tri( iFace );
        }
    }
    else if ( Dim::dimension == ONEFLOW::TWO_D )
    {
        for ( int iFace = 0; iFace < nFaces; ++ iFace )
        {
            this->gFid = iFace;
           // FaceSort * faceSort = this->faceArray[ iFace ];
            //this->SplitLine( faceSort );
            this->SplitLine( iFace );
        }
    }
}

void FaceSearch::SplitQuad2Tri( int faceId )
{
    const IntField & nodeId = this->GetFaceNodeRing( faceId );
    int nNodes = nodeId.size();
    if ( nNodes != 4 ) return;

    LinkField localTriId, localTriFlag;
    this->GetLocalTri( localTriId, localTriFlag );

    LinkField triId;
    this->GetTriId( nodeId, localTriId, triId );

    for ( int iTri = 0; iTri < triId.size(); ++ iTri )
    {
        IntField & tri = triId[ iTri ];
        int fId = this->FindFace( tri );
        if ( fId != ONEFLOW::INVALID_INDEX )
        {
            this->interfaceSplit.childFaceIdsOfParent[ faceId ].push_back( fId );
            this->interfaceSplit.childFaceNodeIds  [ fId ] = localTriId  [ iTri ];
            this->interfaceSplit.childFaceNodeFlags[ fId ] = localTriFlag[ iTri ];
        }
    }
}

void FaceSearch::SplitLine( int faceId )
{
    const IntField & nodeId = this->GetFaceNodeRing( faceId );
    int nNodes = nodeId.size();
    if ( nNodes != 2 ) return;

    LinkField localLineId, localLineFlag;
    LinkField lineId;
    if ( ! this->GetLine( nodeId, localLineId, localLineFlag, lineId ) ) return;

    for ( int iLi = 0; iLi < lineId.size(); ++ iLi )
    {
        IntField & lId = lineId[ iLi ];
        int fId = this->FindFace( lId );
        if ( fId != ONEFLOW::INVALID_INDEX )
        {
            this->interfaceSplit.childFaceIdsOfParent[ faceId ].push_back( fId );
            this->interfaceSplit.childFaceNodeIds  [ fId ] = localLineId  [ iLi ];
            this->interfaceSplit.childFaceNodeFlags[ fId ] = localLineFlag[ iLi ];
        }
    }
}

void BuildQuadDiagonalLocalTris( LinkField & localTriId, LinkField & localTriFlag )
{
    localTriId.clear();
    localTriFlag.clear();

    // Diagonal 0-2
    localTriId.push_back( IntField{ 0, 1, 2 } );
    localTriId.push_back( IntField{ 2, 3, 0 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );

    // Diagonal 1-3
    localTriId.push_back( IntField{ 3, 0, 1 } );
    localTriId.push_back( IntField{ 1, 2, 3 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );
}

void MapLocalTrisToGlobalNodes(
    const IntField & nodeId,
    const LinkField & localTriId,
    LinkField & triId )
{
    triId.clear();
    triId.resize( localTriId.size() );
    for ( std::size_t iTri = 0; iTri < localTriId.size(); ++ iTri )
    {
        for ( int iNode = 0; iNode < 3; ++ iNode )
        {
            const int localCorner = localTriId[ iTri ][ iNode ];
            if ( localCorner < 0 ||
                 static_cast< std::size_t >( localCorner ) >= nodeId.size() )
            {
                throw std::out_of_range(
                    "MapLocalTrisToGlobalNodes: local corner index is out of range" );
            }
            triId[ iTri ].push_back( nodeId[ localCorner ] );
        }
    }
}

void FaceSearch::GetLocalTri( LinkField & localTriId, LinkField & localTriFlag )
{
    BuildQuadDiagonalLocalTris( localTriId, localTriFlag );
}

void FaceSearch::GetTriId( const IntField & nodeId, LinkField & localTriId, LinkField & triId )
{
    MapLocalTrisToGlobalNodes( nodeId, localTriId, triId );
}

bool FaceSearch::GetLine( const IntField & nodeId, LinkField & localLineId, LinkField & localLineFlag, LinkField & lineId )
{
    PointLocator * point_search = this->iFaceLink->point_search.get();
    int p0 = nodeId[ 0 ];
    int p1 = nodeId[ 1 ];
    RealField coor0( 3 ), coor1( 3 ), coor2( 3 );
    point_search->GetPoint( p0, coor0[ 0 ], coor0[ 1 ], coor0[ 2 ] );
    point_search->GetPoint( p1, coor1[ 0 ], coor1[ 1 ], coor1[ 2 ] );
    for ( int m = 0; m < 3; ++ m )
    {
        coor2[ m ] = half * ( coor0[ m ] + coor1[ m ] );
    }
    if ( point_search->FindPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] ) == -1 ) return false;

    const LinkField & globalFaceZones = this->iFaceLink->gI2Zid;
    if ( this->gFid < 0 ||
         static_cast< size_t >( this->gFid ) >= globalFaceZones.size() )
    {
        throw std::out_of_range( "FaceSearch::GetLine: global face ID is out of range" );
    }

    const IntField & zoneIds = globalFaceZones[ this->gFid ];
    if ( zoneIds.size() != 1 )
    {
        throw std::logic_error(
            "FaceSearch::GetLine: a split line must reference exactly one zone" );
    }

    int zoneIndex = zoneIds[ 0 ];
    Grid & grid = this->iFaceLink->GetGrid( zoneIndex );
    int nNodes = grid.nodeMesh->GetNumberOfNodes();
    int pId = nNodes;
    grid.nodeMesh->AddPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] );
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
