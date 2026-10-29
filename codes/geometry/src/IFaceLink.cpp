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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

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


void RegisterLocalInterfaceFaceOnMaps(
    int gIid,
    bool isNew,
    int zid,
    int localFaceId,
    LinkField & l2g,
    LinkField & gI2Zid,
    LinkField & g2l )
{
    if ( zid < 0 || static_cast< size_t >( zid ) >= l2g.size() )
    {
        throw std::out_of_range( "RegisterLocalInterfaceFace: zone index is out of range" );
    }
    if ( localFaceId < 0 || static_cast< size_t >( localFaceId ) >= l2g[ zid ].size() )
    {
        throw std::out_of_range( "RegisterLocalInterfaceFace: local face index is out of range" );
    }
    if ( gI2Zid.size() != g2l.size() )
    {
        throw std::logic_error( "RegisterLocalInterfaceFace: global face mappings are inconsistent" );
    }

    if ( isNew )
    {
        if ( gIid < 0 || static_cast< size_t >( gIid ) != gI2Zid.size() )
        {
            throw std::logic_error( "RegisterLocalInterfaceFace: new global face ID is inconsistent" );
        }

        l2g[ zid ][ localFaceId ] = gIid;
        gI2Zid.push_back( IntField{ zid } );
        g2l.push_back( IntField{ localFaceId } );
    }
    else
    {
        if ( gIid < 0 || static_cast< size_t >( gIid ) >= gI2Zid.size() )
        {
            throw std::out_of_range( "RegisterLocalInterfaceFace: global face ID is out of range" );
        }
        if ( gI2Zid[ gIid ].size() != g2l[ gIid ].size() )
        {
            throw std::logic_error( "RegisterLocalInterfaceFace: global face references are inconsistent" );
        }

        for ( std::size_t slot = 0; slot < gI2Zid[ gIid ].size(); ++ slot )
        {
            if ( gI2Zid[ gIid ][ slot ] == zid &&
                 g2l[ gIid ][ slot ] == localFaceId )
            {
                throw std::logic_error(
                    "RegisterLocalInterfaceFace: local interface face is already registered" );
            }
        }
        if ( gI2Zid[ gIid ].size() >= 2 )
        {
            throw std::logic_error(
                "RegisterLocalInterfaceFace: global face already has two zone participants" );
        }

        l2g[ zid ][ localFaceId ] = gIid;
        gI2Zid[ gIid ].push_back( zid );
        g2l[ gIid ].push_back( localFaceId );
    }
}

IFaceLink::IFaceLink( Grids & gridsIn ) : grids( gridsIn )
{
    const int nZone = GridsSize( gridsIn );
    this->l2g.resize( nZone );

    this->face_search = std::make_unique< FaceSearch >();
    // Create the locator object but do not build its ADT yet. Callers that only
    // need face topology (l2g / face_search) must not require nodeMesh bounds.
    this->point_search = std::make_unique< PointLocator >();
    this->pointLocatorReady = false;
}




void InitNewInterfaceMappingBuffers(
    const LinkField & gI2Zid,
    const LinkField & g2l,
    const LinkField & l2g,
    LinkField & gI2ZidNew,
    LinkField & g2lNew,
    LinkField & l2gNew )
{
    if ( gI2Zid.size() != g2l.size() )
    {
        throw std::logic_error(
            "InitNewInterfaceMappingBuffers: gI2Zid and g2l size mismatch" );
    }

    gI2ZidNew = gI2Zid;
    g2lNew = g2l;

    l2gNew.clear();
    l2gNew.resize( l2g.size() );
}

void CommitInterfaceMappingBuffers(
    LinkField & gI2Zid,
    LinkField & g2l,
    LinkField & l2g,
    LinkField & gI2ZidNew,
    LinkField & g2lNew,
    LinkField & l2gNew,
    const char * operation )
{
    ValidateGlobalFaceMappings( gI2ZidNew, g2lNew, operation );

    gI2Zid.swap( gI2ZidNew );
    g2l.swap( g2lNew );
    l2g.swap( l2gNew );
}

int CountInterfaceFacesMissingPartner(
    const IntField & localToGlobal,
    const LinkField & zoneIdsPerGlobalFace,
    const LinkField & localFaceIdsPerGlobalFace )
{
    int missing = 0;
    for ( std::size_t i = 0; i < localToGlobal.size(); ++ i )
    {
        const int gIFace = localToGlobal[ i ];
        ValidateGlobalFaceMapping(
            gIFace,
            zoneIdsPerGlobalFace,
            localFaceIdsPerGlobalFace,
            "CountInterfaceFacesMissingPartner" );
        const int nIZone = static_cast< int >( zoneIdsPerGlobalFace[ gIFace ].size() );
        // Legacy: only count faces with 0 or 1 zone slots as "missing".
        // Over-subscribed faces (nIZone > 2) are ignored in the counter.
        if ( nIZone < 2 )
        {
            ++ missing;
        }
    }
    return missing;
}

bool FindPeriodicInterfacePartner(
    int selfZoneId,
    int selfLocalFaceId,
    const IntField & zoneIds,
    const IntField & localFaceIds,
    int & outRemoteZoneId,
    int & outRemoteLocalFaceId )
{
    if ( zoneIds.size() != localFaceIds.size() )
    {
        throw std::logic_error(
            "FindPeriodicInterfacePartner: zone and local-face rows are inconsistent" );
    }

    // Periodic matching searches the geometrically transformed face, whose
    // participant row may not contain the current face. Skip self when present,
    // but do not require self membership as FindRemoteInterfacePartner does.
    for ( std::size_t i = 0; i < zoneIds.size(); ++ i )
    {
        if ( zoneIds[ i ] == selfZoneId && localFaceIds[ i ] == selfLocalFaceId )
        {
            continue;
        }
        outRemoteZoneId = zoneIds[ i ];
        outRemoteLocalFaceId = localFaceIds[ i ];
        return true;
    }
    return false;
}

bool FindRemoteInterfacePartner(
    int globalFaceId,
    int selfZoneId,
    int selfLocalFaceId,
    const LinkField & zoneIdsPerGlobalFace,
    const LinkField & localFaceIdsPerGlobalFace,
    int & outRemoteZoneId,
    int & outRemoteLocalFaceId )
{
    ValidateGlobalFaceMapping(
        globalFaceId,
        zoneIdsPerGlobalFace,
        localFaceIdsPerGlobalFace,
        "FindRemoteInterfacePartner" );

    const IntField & zoneIds = zoneIdsPerGlobalFace[ globalFaceId ];
    const IntField & localIds = localFaceIdsPerGlobalFace[ globalFaceId ];

    // A remote participant is meaningful only if the caller is registered
    // on this global face. Otherwise the tables do not describe this lookup.
    bool foundSelf = false;
    for ( std::size_t i = 0; i < zoneIds.size(); ++ i )
    {
        if ( zoneIds[ i ] == selfZoneId && localIds[ i ] == selfLocalFaceId )
        {
            foundSelf = true;
            break;
        }
    }
    if ( ! foundSelf )
    {
        return false;
    }

    for ( std::size_t i = 0; i < zoneIds.size(); ++ i )
    {
        if ( zoneIds[ i ] != selfZoneId || localIds[ i ] != selfLocalFaceId )
        {
            outRemoteZoneId = zoneIds[ i ];
            outRemoteLocalFaceId = localIds[ i ];
            return true;
        }
    }
    return false;
}

void IFaceLink::EnsurePointLocatorReady()
{
    if ( this->pointLocatorReady )
    {
        return;
    }
    if ( ! this->point_search )
    {
        this->point_search = std::make_unique< PointLocator >();
    }
    this->point_search->Initialize( this->grids );
    this->pointLocatorReady = true;
}

PointLocator & IFaceLink::GetPointLocator()
{
    if ( ! this->point_search )
    {
        throw std::logic_error(
            "IFaceLink::GetPointLocator: point locator is not initialized" );
    }
    return *this->point_search;
}

const PointLocator & IFaceLink::GetPointLocator() const
{
    if ( ! this->point_search )
    {
        throw std::logic_error(
            "IFaceLink::GetPointLocator: point locator is not initialized" );
    }
    return *this->point_search;
}

void IFaceLink::ReplacePointLocator( std::unique_ptr< PointLocator > locator, bool markReady )
{
    if ( ! locator )
    {
        throw std::invalid_argument(
            "IFaceLink::ReplacePointLocator: locator must not be null" );
    }
    this->point_search = std::move( locator );
    this->pointLocatorReady = markReady;
}

FaceSearch & IFaceLink::GetFaceSearch()
{
    if ( ! this->face_search )
    {
        throw std::logic_error(
            "IFaceLink::GetFaceSearch: face search is not initialized" );
    }
    return *this->face_search;
}

const FaceSearch & IFaceLink::GetFaceSearch() const
{
    if ( ! this->face_search )
    {
        throw std::logic_error(
            "IFaceLink::GetFaceSearch: face search is not initialized" );
    }
    return *this->face_search;
}

InterfaceFaceSplitTable & IFaceLink::MutableInterfaceSplit()
{
    return this->GetFaceSearch().interfaceSplit;
}

const InterfaceFaceSplitTable & IFaceLink::GetInterfaceSplit() const
{
    return this->GetFaceSearch().interfaceSplit;
}



IFaceLink::~IFaceLink() = default;

Grid & IFaceLink::GetGrid( int zoneIndex )
{
    if ( zoneIndex < 0 || static_cast< size_t >( zoneIndex ) >= this->grids.size() )
    {
        throw std::out_of_range( "IFaceLink::GetGrid: zone index is out of range" );
    }
    return GridAt( this->grids, static_cast< size_t >( zoneIndex ) );
}

const Grid & IFaceLink::GetGrid( int zoneIndex ) const
{
    if ( zoneIndex < 0 || static_cast< size_t >( zoneIndex ) >= this->grids.size() )
    {
        throw std::out_of_range( "IFaceLink::GetGrid: zone index is out of range" );
    }
    return GridAt( this->grids, static_cast< size_t >( zoneIndex ) );
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
    this->AddFace( faceNode );

    auto [gIid, isNew] = this->faceLookup.FindOrAdd( faceNode );
    RegisterLocalInterfaceFaceOnMaps(
        gIid, isNew, zid, lCount, this->l2g, this->gI2Zid, this->g2l );
}

void IFaceLink::ReconstructInterFace()
{
    this->face_search->CalcNewFaceId( *this );
}

void IFaceLink::InitNewLgMapping()
{
    InitNewInterfaceMappingBuffers(
        this->gI2Zid,
        this->g2l,
        this->l2g,
        this->gI2ZidNew,
        this->g2lNew,
        this->l2gNew );
}

void IFaceLink::UpdateLgMapping()
{
    CommitInterfaceMappingBuffers(
        this->gI2Zid,
        this->g2l,
        this->l2g,
        this->gI2ZidNew,
        this->g2lNew,
        this->l2gNew,
        "IFaceLink::UpdateLgMapping" );
}

void IFaceLink::MatchInterfaceTopology( Grid & grid )
{
    ValidateGridIndex( grid, "IFaceLink::MatchInterfaceTopology" );
    InterFace * interFace = grid.interFace.get();
    if ( ! interFace ) return;

    const IntField & localToGlobal = this->l2g[ grid.id ];
    const int missingPartnerCount = CountInterfaceFacesMissingPartner(
        localToGlobal, this->gI2Zid, this->g2l );

    const int nIFaces = static_cast< int >( localToGlobal.size() );
    for ( int iIFace = 0; iIFace < nIFaces; ++ iIFace )
    {
        const int gIFace = localToGlobal[ iIFace ];
        // Recompute both fields on every pass. A face that was previously
        // paired may become unpaired after interface-map rebuilding.
        interFace->zoneId[ iIFace ] = -1;
        interFace->localInterfaceId[ iIFace ] = -1;

        int remoteZone = -1;
        int remoteLocal = -1;
        if ( FindRemoteInterfacePartner(
                 gIFace,
                 grid.id,
                 iIFace,
                 this->gI2Zid,
                 this->g2l,
                 remoteZone,
                 remoteLocal ) )
        {
            interFace->zoneId[ iIFace ] = remoteZone;
            interFace->localInterfaceId[ iIFace ] = remoteLocal;
        }
    }
    std::cout << " Interface faces missing a partner = "
              << missingPartnerCount << "\n";
}


void IFaceLink::MatchPeriodicInterface( Grid & grid )
{
    ValidateGridIndex( grid, "IFaceLink::MatchPeriodicInterface" );
    InterFace * interFace = grid.interFace.get();
    if ( ! interFace ) return;

    this->EnsurePointLocatorReady();

    int nIFaces = this->l2g[ grid.id ].size();

    for ( int iIFace = 0; iIFace < nIFaces; ++ iIFace )
    {
        int gIFace = this->l2g[ grid.id ][ iIFace ];
        ValidateGlobalFaceMapping(
            gIFace, this->gI2Zid, this->g2l, "IFaceLink::MatchPeriodicInterface" );
        int nIZone = static_cast< int >( this->gI2Zid[ gIFace ].size() );

        // A two-zone face already has a topological partner. Other faces
        // must be rematched geometrically, so clear any result from an
        // earlier pass before any lookup that may fail.
        if ( nIZone == 2 ) continue;
        interFace->zoneId[ iIFace ] = -1;
        interFace->localInterfaceId[ iIFace ] = -1;

        // GetFaceNodeRing bounds-checks against stored face count.
        const IntField & nodeId = this->face_search->GetFaceNodeRing( gIFace );

        RealField xList, yList, zList;
        this->GetPointLocator().GetFaceCoorList( nodeId, xList, yList, zList );

        RealField xxList, yyList, zzList;
        f2fmap.FindFace( xList, yList, zList, xxList, yyList, zzList );

        int nNodes = xxList.size();
        IntField faceNode_period;

        for ( int i = 0; i < nNodes; ++ i )
        {
            Real xm = xxList[ i ];
            Real ym = yyList[ i ];
            Real zm = zzList[ i ];
            int id = this->GetPointLocator().FindPoint( xm, ym, zm );
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

        // Ensure periodic partner is a known stored face before zone lookup.
        (void) this->face_search->GetFaceNodeRing( faceId_period );

        const IntField & periodicZones = this->gI2Zid[ faceId_period ];
        const IntField & periodicLocalIds = this->g2l[ faceId_period ];
        if ( periodicZones.size() != periodicLocalIds.size() )
        {
            throw std::logic_error(
                "IFaceLink::MatchPeriodicInterface: periodic face references are inconsistent" );
        }
        if ( periodicZones.empty() ) continue;

        int nZid_period = -1;
        int lId_period = -1;
        if ( ! FindPeriodicInterfacePartner(
                 grid.id, iIFace, periodicZones, periodicLocalIds,
                 nZid_period, lId_period ) )
        {
            continue;
        }

        interFace->zoneId[ iIFace ] = nZid_period;
        interFace->localInterfaceId[ iIFace ] = lId_period;
    }
}

void GetFaceCoorList( const IntField & faceNode, RealField & xList, RealField & yList, RealField & zList, const NodeMesh & nodeMesh )
{
    const size_t nPoint = faceNode.size();
    if ( xList.size() < nPoint || yList.size() < nPoint || zList.size() < nPoint )
    {
        throw std::invalid_argument( "GetFaceCoorList: coordinate buffers are too small" );
    }

    const size_t nNodes = nodeMesh.xN.size();
    if ( nodeMesh.yN.size() != nNodes || nodeMesh.zN.size() != nNodes )
    {
        throw std::logic_error( "GetFaceCoorList: node coordinate arrays have different sizes" );
    }

    for ( size_t iNode = 0; iNode < nPoint; ++ iNode )
    {
        const int gN = faceNode[ iNode ];
        if ( gN < 0 || static_cast< size_t >( gN ) >= nNodes )
        {
            throw std::out_of_range( "GetFaceCoorList: node index is out of range" );
        }

        xList[ iNode ] = nodeMesh.xN[ gN ];
        yList[ iNode ] = nodeMesh.yN[ gN ];
        zList[ iNode ] = nodeMesh.zN[ gN ];
    }
}

void GetCoorIdList( IFaceLink & iFaceLink, RealField & xList, RealField & yList, RealField & zList, int nPoint, IntField & pointId )
{
    if ( nPoint < 0 )
    {
        throw std::invalid_argument( "GetCoorIdList: point count must not be negative" );
    }

    const size_t count = static_cast< size_t >( nPoint );
    if ( xList.size() < count || yList.size() < count ||
         zList.size() < count || pointId.size() < count )
    {
        throw std::invalid_argument( "GetCoorIdList: input or output buffers are too small" );
    }
    iFaceLink.EnsurePointLocatorReady();

    for ( size_t iNode = 0; iNode < count; ++ iNode )
    {
        const Real xm = xList[ iNode ];
        const Real ym = yList[ iNode ];
        const Real zm = zList[ iNode ];

        pointId[ iNode ] = iFaceLink.GetPointLocator().AddPoint( xm, ym, zm );
    }
}



void IFaceLink::EnsureInterfaceZoneCount( std::size_t nZones )
{
    if ( this->l2g.size() < nZones )
    {
        this->l2g.resize( nZones );
    }
    if ( this->l2gNew.size() < nZones )
    {
        this->l2gNew.resize( nZones );
    }
}

void IFaceLink::EnsureLocalInterfaceFaceCount( int zoneId, std::size_t localFaceCount )
{
    if ( zoneId < 0 || static_cast< size_t >( zoneId ) >= this->l2g.size() )
    {
        throw std::out_of_range(
            "IFaceLink::EnsureLocalInterfaceFaceCount: zone id is out of range" );
    }
    if ( this->l2g[ zoneId ].size() < localFaceCount )
    {
        this->l2g[ zoneId ].resize( localFaceCount, -1 );
    }
}

void IFaceLink::RegisterLocalInterfaceFace( int gIid, bool isNew, int zid, int localFaceId )
{
    RegisterLocalInterfaceFaceOnMaps(
        gIid, isNew, zid, localFaceId, this->l2g, this->gI2Zid, this->g2l );
}

void IFaceLink::EnsureSplitBookkeepingMatchesZones()
{
    // nChild is parallel to zones in l2g (not to global faces).
    this->EnsureInterfaceZoneCount( this->l2g.size() );
    if ( this->nChild.size() < this->l2g.size() )
    {
        this->nChild.resize( this->l2g.size() );
    }
}

IntField & IFaceLink::MutableL2gNewForZone( int zoneId )
{
    if ( zoneId < 0 || static_cast< size_t >( zoneId ) >= this->l2gNew.size() )
    {
        throw std::out_of_range(
            "IFaceLink::MutableL2gNewForZone: zone id is out of range for l2gNew" );
    }
    return this->l2gNew[ zoneId ];
}

IntField & IFaceLink::MutableNChildForZone( int zoneId )
{
    if ( zoneId < 0 || static_cast< size_t >( zoneId ) >= this->nChild.size() )
    {
        throw std::out_of_range(
            "IFaceLink::MutableNChildForZone: zone id is out of range for nChild" );
    }
    return this->nChild[ zoneId ];
}

const IntField & IFaceLink::GetL2gForZone( int zoneId ) const
{
    if ( zoneId < 0 || static_cast< size_t >( zoneId ) >= this->l2g.size() )
    {
        throw std::out_of_range(
            "IFaceLink::GetL2gForZone: zone id is out of range for l2g" );
    }
    return this->l2g[ zoneId ];
}

EndNameSpace
