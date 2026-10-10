/*---------------------------------------------------------------------------*\\
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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\\*---------------------------------------------------------------------------*/

#include "IFaceLink.h"
#include "Constant.h"
#include "InterFace.h"
#include "Grid.h"
#include "PointLocator.h"
#include "FaceSearch.h"
#include "CgnsPeriod.h"
#include "NodeMesh.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{
void ValidateGlobalFaceMapping(
    int globalFaceId, const LinkField & zoneIds, const LinkField & localFaceIds,
    const char * operation )
{
    if ( globalFaceId < 0 ||
         static_cast< size_t >( globalFaceId ) >= zoneIds.size() ||
         static_cast< size_t >( globalFaceId ) >= localFaceIds.size() )
    {
        throw std::logic_error( std::string( operation ) + ": global face mapping is out of range" );
    }

    if ( zoneIds[ globalFaceId ].size() != localFaceIds[ globalFaceId ].size() )
    {
        throw std::logic_error( std::string( operation ) + ": global face references are inconsistent" );
    }
}

void ValidateGlobalFaceMappings(
    const LinkField & zoneIds, const LinkField & localFaceIds,
    const char * operation )
{
    if ( zoneIds.size() != localFaceIds.size() )
    {
        throw std::logic_error( std::string( operation ) + ": global face mapping tables have different sizes" );
    }

    for ( size_t globalFaceId = 0; globalFaceId < zoneIds.size(); ++ globalFaceId )
    {
        ValidateGlobalFaceMapping(
            static_cast< int >( globalFaceId ), zoneIds, localFaceIds, operation );
    }
}
}

IFaceLink::IFaceLink( Grids & gridsIn ) : grids( gridsIn )
{
    const int nZone = GridsSize( gridsIn );
    this->l2g.resize( nZone );

    this->face_search = std::make_unique< FaceSearch >();
    this->point_search = std::make_unique< PointLocator >();
    this->point_search->Initialize( gridsIn );
}

IFaceLink::~IFaceLink() = default;

Grid & IFaceLink::GetGrid( int zoneIndex )
{
    return GridAt( this->grids, zoneIndex );
}

void IFaceLink::ValidateGridIndex( const Grid & grid, const char * operation ) const
{
    const int zid = grid.id;
    if ( zid < 0 || static_cast< size_t >( zid ) >= this->l2g.size() )
    {
        throw std::out_of_range( std::string( operation ) + ": grid zone index is out of range" );
    }
    if ( & GridAt( this->grids, static_cast< size_t >( zid ) ) != & grid )
    {
        throw std::invalid_argument(
            std::string( operation ) + ": grid does not belong to the linked collection" );
    }
}

void IFaceLink::Init( Grid & grid )
{
    ValidateGridIndex( grid, "IFaceLink::Init" );
    const int zid = grid.id;
    if ( ! grid.interFace )
    {
        throw std::logic_error( "IFaceLink::Init: grid interface data is not initialized" );
    }

    const int nIFaces = grid.interFace->nIFaces;
    if ( nIFaces < 0 )
    {
        throw std::invalid_argument( "IFaceLink::Init: interface face count must not be negative" );
    }

    this->l2g[ zid ].resize( static_cast< size_t >( nIFaces ) );
}

void IFaceLink::AddFace( const IntField & facePointIndexes )
{
    this->face_search->AddFace( facePointIndexes );
}

void IFaceLink::CreateLink( IntField & faceNode, int zid, int lCount )
{
    if ( zid < 0 || static_cast< size_t >( zid ) >= this->l2g.size() )
    {
        throw std::out_of_range( "IFaceLink::CreateLink: zone index is out of range" );
    }
    if ( lCount < 0 || static_cast< size_t >( lCount ) >= this->l2g[ zid ].size() )
    {
        throw std::out_of_range( "IFaceLink::CreateLink: local face index is out of range" );
    }
    if ( this->gI2Zid.size() != this->g2l.size() )
    {
        throw std::logic_error( "IFaceLink::CreateLink: global face mappings are inconsistent" );
    }

    this->AddFace( faceNode );

    auto [gIid, isNew] = this->faceLookup.FindOrAdd( faceNode );

    if ( isNew )
    {
        if ( gIid < 0 || static_cast< size_t >( gIid ) != this->gI2Zid.size() )
        {
            throw std::logic_error( "IFaceLink::CreateLink: new global face ID is inconsistent" );
        }

        this->l2g[ zid ][ lCount ] = gIid;

        IntField zids;
        IntField lIid;
        zids.push_back( zid );
        lIid.push_back( lCount );
        this->gI2Zid.push_back( std::move( zids ) );
        this->g2l.push_back( std::move( lIid ) );
    }
    else
    {
        if ( gIid < 0 || static_cast< size_t >( gIid ) >= this->gI2Zid.size() )
        {
            throw std::logic_error( "IFaceLink::CreateLink: global face ID is out of range" );
        }
        if ( this->gI2Zid[ gIid ].size() != this->g2l[ gIid ].size() )
        {
            throw std::logic_error( "IFaceLink::CreateLink: global face references are inconsistent" );
        }

        this->l2g[ zid ][ lCount ] = gIid;
        this->gI2Zid[ gIid ].push_back( zid );
        this->g2l[ gIid ].push_back( lCount );
    }
}

void IFaceLink::ReconstructInterFace()
{
    this->face_search->CalcNewFaceId( *this );
}

void IFaceLink::InitNewLgMapping()
{
    this->gI2ZidNew = this->gI2Zid;
    this->g2lNew = this->g2l;
}

void IFaceLink::UpdateLgMapping()
{
    ValidateGlobalFaceMappings(
        this->gI2ZidNew, this->g2lNew, "IFaceLink::UpdateLgMapping" );

    this->gI2Zid.swap( this->gI2ZidNew );
    this->g2l.swap( this->g2lNew );
}

void IFaceLink::MatchInterfaceTopology( Grid & grid )
{
    ValidateGridIndex( grid, "IFaceLink::MatchInterfaceTopology" );
    InterFace * interFace = grid.interFace.get();
    if ( ! interFace ) return;

    int missingPeriodicPartnerCount = 0;

    int nIFaces = this->l2g[ grid.id ].size();

    for ( int iIFace = 0; iIFace < nIFaces; ++ iIFace )
    {
        int gIFace = this->l2g[ grid.id ][ iIFace ];
        ValidateGlobalFaceMapping(
            gIFace, this->gI2Zid, this->g2l, "IFaceLink::MatchInterfaceTopology" );
        int nIZone = static_cast< int >( this->gI2Zid[ gIFace ].size() );

        if ( nIZone != 2 )
        {
            if ( nIZone > 2 )
            {
            }
            else
            {
                ++ missingPeriodicPartnerCount;
            }
        }

        for ( int iIZone = 0; iIZone < nIZone; ++ iIZone )
        {
            int nZid = this->gI2Zid[ gIFace ][ iIZone ];
            int lId  = this->g2l[ gIFace ][ iIZone ];
            if ( ( nZid != grid.id ) ||
                 ( lId  != iIFace   ) )
            {
                interFace->zoneId[ iIFace ] = nZid;
                interFace->localInterfaceId[ iIFace ] = lId;
                break;
            }
        }
    }
    std::cout << " Periodic boundary faces missing a partner = "
              << missingPeriodicPartnerCount << "\n";
}

void IFaceLink::MatchPeriodicInterface( Grid & grid )
{
    ValidateGridIndex( grid, "IFaceLink::MatchPeriodicInterface" );
    InterFace * interFace = grid.interFace.get();
    if ( ! interFace ) return;

    int nIFaces = this->l2g[ grid.id ].size();

    for ( int iIFace = 0; iIFace < nIFaces; ++ iIFace )
    {
        int gIFace = this->l2g[ grid.id ][ iIFace ];
        ValidateGlobalFaceMapping(
            gIFace, this->gI2Zid, this->g2l, "IFaceLink::MatchPeriodicInterface" );
        if ( static_cast< size_t >( gIFace ) >= this->face_search->faceArray.size() )
        {
            throw std::logic_error( "IFaceLink::MatchPeriodicInterface: global face ID is out of range" );
        }
        int nIZone = static_cast< int >( this->gI2Zid[ gIFace ].size() );

        if ( nIZone == 2 ) continue;

        const IntField & nodeId = this->face_search->faceArray[ gIFace ];

        RealField xList, yList, zList;
        this->point_search->GetFaceCoorList( nodeId, xList, yList, zList );

        RealField xxList, yyList, zzList;
        f2fmap.FindFace( xList, yList, zList, xxList, yyList, zzList );

        int nNodes = xxList.size();
        IntField faceNode_period;

        for ( int i = 0; i < nNodes; ++ i )
        {
            Real xm = xxList[ i ];
            Real ym = yyList[ i ];
            Real zm = zzList[ i ];
            int id = this->point_search->FindPoint( xm, ym, zm );
            if ( id == INVALID_INDEX )
            {
                faceNode_period.clear();
                break;
            }
            faceNode_period.push_back( id );
        }

        if ( faceNode_period.size() != nNodes ) continue;

        int faceId_period = this->face_search->FindFace( faceNode_period );
        if ( faceId_period == INVALID_INDEX ) continue;
        ValidateGlobalFaceMapping(
            faceId_period, this->gI2Zid, this->g2l,
            "IFaceLink::MatchPeriodicInterface periodic partner" );

        if ( static_cast< size_t >( faceId_period ) >= this->face_search->faceArray.size() )
        {
            throw std::logic_error(
                "IFaceLink::MatchPeriodicInterface: periodic face ID is out of range" );
        }

        const IntField & periodicZones = this->gI2Zid[ faceId_period ];
        const IntField & periodicLocalIds = this->g2l[ faceId_period ];
        if ( periodicZones.size() != periodicLocalIds.size() )
        {
            throw std::logic_error(
                "IFaceLink::MatchPeriodicInterface: periodic face references are inconsistent" );
        }
        if ( periodicZones.empty() ) continue;

        int nZid_period = periodicZones[ 0 ];
        int lId_period  = periodicLocalIds[ 0 ];

        interFace->zoneId[ iIFace ] = nZid_period;
        interFace->localInterfaceId[ iIFace ] = lId_period;
    }
}

void GetFaceCoorList( const IntField & faceNode, RealField & xList, RealField & yList, RealField & zList, const NodeMesh & nodeMesh )
{
    int nPoint = faceNode.size();
    for ( int iNode = 0; iNode < nPoint; ++ iNode )
    {
        int gN = faceNode[ iNode ];
        xList[ iNode ] = nodeMesh.xN[ gN ];
        yList[ iNode ] = nodeMesh.yN[ gN ];
        zList[ iNode ] = nodeMesh.zN[ gN ];
    }
}

void GetCoorIdList( IFaceLink & iFaceLink, RealField & xList, RealField & yList, RealField & zList, int nPoint, IntField & pointId )
{
    for ( int iNode = 0; iNode < nPoint; ++ iNode )
    {
        Real xm = xList[ iNode ];
        Real ym = yList[ iNode ];
        Real zm = zList[ iNode ];

        pointId[ iNode ] = iFaceLink.point_search->AddPoint( xm, ym, zm );
    }
}

EndNameSpace
