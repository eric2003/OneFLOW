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

#include "ScalarIFace.h"
#include "MetisGrid.h"
#include "DataStorage.h"
#include "DataBaseIO.h"
#include "DataBook.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <utility>


BeginNameSpace( ONEFLOW )

ScalarIFaceIJ::ScalarIFaceIJ()
{
    ;
}

ScalarIFaceIJ::~ScalarIFaceIJ()
{
}

void ScalarIFaceIJ::WriteInterfaceTopology( DataBook * databook )
{
    if ( this->ifaces.size() != this->recv_ifaces.size() )
    {
        throw std::logic_error( "ScalarIFaceIJ::WriteInterfaceTopology: interface arrays have inconsistent sizes" );
    }

    int nIFaces = static_cast< int >( this->ifaces.size() );
    ONEFLOW::HXWrite( databook, this->zonej );
    ONEFLOW::HXWrite( databook, nIFaces );
    ONEFLOW::HXWrite( databook, this->ifaces );
    ONEFLOW::HXWrite( databook, this->recv_ifaces );
}

void ScalarIFaceIJ::ReadInterfaceTopology( DataBook * databook )
{
    int zonej = -1;
    int nIFaces = -1;
    ONEFLOW::HXRead( databook, zonej );
    ONEFLOW::HXRead( databook, nIFaces );
    if ( nIFaces < 0 )
    {
        throw std::runtime_error( "ScalarIFaceIJ::ReadInterfaceTopology: interface count must be non-negative" );
    }

    std::vector< int > ifaces( nIFaces );
    std::vector< int > recvIfaces( nIFaces );
    ONEFLOW::HXRead( databook, ifaces );
    ONEFLOW::HXRead( databook, recvIfaces );

    // Commit the decoded neighbor topology only after all fields have been read.
    this->zonej = zonej;
    this->ifaces = std::move( ifaces );
    this->recv_ifaces = std::move( recvIfaces );
}

ScalarIFace::ScalarIFace()
    : dataSend( std::make_unique< DataStorage >() ),
      dataRecv( std::make_unique< DataStorage >() )
{
}

ScalarIFace::~ScalarIFace() = default;

void ScalarIFace::AddInterface( int global_interface_id, int neighbor_zoneid, int neighbor_cellid )
{
    if ( global_interface_id < 0 || neighbor_zoneid < 0 || neighbor_cellid < 0 )
    {
        throw std::invalid_argument( "ScalarIFace::AddInterface: interface and neighbor IDs must be non-negative" );
    }

    const size_t nInterfaces = this->iglobalfaces.size();
    if ( this->zones.size() != nInterfaces || this->cells.size() != nInterfaces ||
         this->global_to_local_interfaces.size() != nInterfaces ||
         this->local_to_global_interfaces.size() != nInterfaces )
    {
        throw std::logic_error( "ScalarIFace::AddInterface: existing interface mappings are inconsistent" );
    }
    if ( this->global_to_local_interfaces.find( global_interface_id ) != this->global_to_local_interfaces.end() )
    {
        throw std::invalid_argument( "ScalarIFace::AddInterface: duplicate global interface ID" );
    }

    const int localInterfaceId = static_cast< int >( nInterfaces );

    // Allocate vector capacity before changing the logical interface mapping.
    this->iglobalfaces.reserve( nInterfaces + 1 );
    this->zones.reserve( nInterfaces + 1 );
    this->cells.reserve( nInterfaces + 1 );

    const auto globalEntry = this->global_to_local_interfaces.emplace( global_interface_id, localInterfaceId );
    if ( ! globalEntry.second )
    {
        throw std::invalid_argument( "ScalarIFace::AddInterface: duplicate global interface ID" );
    }

    try
    {
        const auto localEntry = this->local_to_global_interfaces.emplace( localInterfaceId, global_interface_id );
        if ( ! localEntry.second )
        {
            throw std::logic_error( "ScalarIFace::AddInterface: local interface ID is already mapped" );
        }
    }
    catch ( ... )
    {
        this->global_to_local_interfaces.erase( globalEntry.first );
        throw;
    }

    // Integer appends cannot allocate after the reserves above, so all five
    // representations are committed together.
    this->iglobalfaces.push_back( global_interface_id );
    this->zones.push_back( neighbor_zoneid );
    this->cells.push_back( neighbor_cellid );
}

int ScalarIFace::GetLocalInterfaceId( int global_interface_id )
{
    const auto iter = this->global_to_local_interfaces.find( global_interface_id );
    if ( iter == this->global_to_local_interfaces.end() )
    {
        throw std::runtime_error( "ScalarIFace::GetLocalInterfaceId: global interface id was not found" );
    }
    return iter->second;
}

int ScalarIFace::GetNIFaces()
{
    return zones.size();
}

int ScalarIFace::FindINeibor( int iZone )
{
    int nNeis = data.size();
    for( int iNei = 0; iNei < nNeis; ++ iNei)
    {
        ScalarIFaceIJ & iFaceIJ = data[ iNei ];
        if ( iZone == iFaceIJ.zonej )
        {
            return iNei;
        }
    }
    return -1;
}

void ScalarIFace::CalcLocalInterfaceId( int iZone, std::vector<int> & globalfaces, std::vector<int> & localfaces )
{
    std::vector< int > reconstructedLocalFaces;
    reconstructedLocalFaces.reserve( globalfaces.size() );
    for ( int i = 0; i < globalfaces.size(); ++ i )
    {
        const int gid = globalfaces[ i ];
        const auto iter = this->global_to_local_interfaces.find( gid );
        if ( iter == this->global_to_local_interfaces.end() )
        {
            throw std::runtime_error( "ScalarIFace::CalcLocalInterfaceId: global interface id was not found" );
        }
        reconstructedLocalFaces.push_back( iter->second );
    }

    // The neighbor of iZone must have a reciprocal entry in this interface list.
    const int jNei = FindINeibor( iZone );
    if ( jNei < 0 )
    {
        throw std::runtime_error( "ScalarIFace::CalcLocalInterfaceId: reciprocal neighbor zone was not found" );
    }

    // Replace derived mappings instead of appending duplicate IDs on repeated reconstruction.
    localfaces = reconstructedLocalFaces;
    this->data[ jNei ].recv_ifaces = std::move( reconstructedLocalFaces );
}

void ScalarIFace::DumpInterfaceMap()
{
    std::cout << " global_to_local_interfaces std::map \n";
    this->DumpMap( this->global_to_local_interfaces );
    std::cout << "\n";
    std::cout << " local_to_global_interfaces std::map \n";
    this->DumpMap( this->local_to_global_interfaces );
    std::cout << "\n";
}

void ScalarIFace::DumpMap( std::map<int,int> & mapin )
{
    for ( std::map<int, int>::iterator iter = mapin.begin(); iter != mapin.end(); ++ iter )
    {
        std::cout << iter->first << " " << iter->second << "\n";
    }
    std::cout << "\n";
}

void ScalarIFace::ReconstructNeighbor()
{
    const size_t nInterfaces = zones.size();
    if ( cells.size() != nInterfaces || iglobalfaces.size() != nInterfaces )
    {
        throw std::runtime_error( "ScalarIFace::ReconstructNeighbor: interface mapping arrays have inconsistent sizes" );
    }

    std::set<int> neighborZones;
    for ( const int neighborZone : zones )
    {
        if ( neighborZone < 0 )
        {
            throw std::runtime_error( "ScalarIFace::ReconstructNeighbor: neighbor zone id must be non-negative" );
        }
        neighborZones.insert( neighborZone );
    }

    std::vector< ScalarIFaceIJ > reconstructed;
    reconstructed.reserve( neighborZones.size() );
    for ( const int neighborZone : neighborZones )
    {
        ScalarIFaceIJ interfaceData;
        interfaceData.zonej = neighborZone;

        for ( size_t iInterface = 0; iInterface < nInterfaces; ++ iInterface )
        {
            if ( zones[ iInterface ] == neighborZone )
            {
                interfaceData.cells.push_back( cells[ iInterface ] );
                interfaceData.iglobalfaces.push_back( iglobalfaces[ iInterface ] );
                interfaceData.ifaces.push_back( static_cast< int >( iInterface ) );
            }
        }
        reconstructed.push_back( std::move( interfaceData ) );
    }

    // Replace derived neighbor data so repeated reconstruction cannot append duplicates.
    data = std::move( reconstructed );
}

void ScalarIFace::WriteInterfaceTopology( DataBook * databook )
{
    const size_t nInterfaces = this->zones.size();
    if ( nInterfaces > 0 &&
         ( this->target_interfaces.size() != nInterfaces ||
           this->interface_to_bcface.size() != nInterfaces ) )
    {
        throw std::logic_error( "ScalarIFace::WriteInterfaceTopology: interface arrays have inconsistent sizes" );
    }
    for ( const ScalarIFaceIJ & neighborData : this->data )
    {
        if ( neighborData.ifaces.size() != neighborData.recv_ifaces.size() )
        {
            throw std::logic_error( "ScalarIFace::WriteInterfaceTopology: neighbor interface arrays have inconsistent sizes" );
        }
    }

    int nIFaces = static_cast< int >( nInterfaces );
    ONEFLOW::HXWrite( databook, nIFaces );
    if ( nIFaces > 0 )
    {
    	ONEFLOW::HXWrite( databook, this->zones               );
    	ONEFLOW::HXWrite( databook, this->target_interfaces   );
    	ONEFLOW::HXWrite( databook, this->interface_to_bcface );

        int nNeis = data.size();
        ONEFLOW::HXWrite( databook, nNeis );
        for ( int iNei = 0; iNei < nNeis; ++ iNei )
        {
            ScalarIFaceIJ & iFaceIJ = data[ iNei ];
            iFaceIJ.WriteInterfaceTopology( databook );
        }
    }
}

void ScalarIFace::ReadInterfaceTopology( DataBook * databook )
{
    int nIFaces = -1;
    ONEFLOW::HXRead( databook, nIFaces );
    if ( nIFaces < 0 )
    {
        throw std::runtime_error( "ScalarIFace::ReadInterfaceTopology: interface count must be non-negative" );
    }

    std::cout << " nIFaces = " << nIFaces << std::endl;

    std::vector< int > zones;
    std::vector< int > targetInterfaces;
    std::vector< int > interfaceToBcface;
    std::vector< ScalarIFaceIJ > interfaceData;

    if ( nIFaces > 0 )
    {
        zones.resize( nIFaces );
        targetInterfaces.resize( nIFaces );
        interfaceToBcface.resize( nIFaces );

        ONEFLOW::HXRead( databook, zones );
        ONEFLOW::HXRead( databook, targetInterfaces );
        ONEFLOW::HXRead( databook, interfaceToBcface );

        int nNeis = -1;
        ONEFLOW::HXRead( databook, nNeis );
        if ( nNeis < 0 )
        {
            throw std::runtime_error( "ScalarIFace::ReadInterfaceTopology: neighbor count must be non-negative" );
        }

        interfaceData.resize( nNeis );
        for ( int iNei = 0; iNei < nNeis; ++ iNei )
        {
            interfaceData[ iNei ].ReadInterfaceTopology( databook );
        }
    }

    // Replace serialized interface state only after the complete read succeeds.
    this->zones = std::move( zones );
    this->target_interfaces = std::move( targetInterfaces );
    this->interface_to_bcface = std::move( interfaceToBcface );
    this->data = std::move( interfaceData );

    // These mappings are not serialized; retaining them would associate the new
    // topology with interface IDs from the previously loaded mesh.
    this->iglobalfaces.clear();
    this->cells.clear();
    this->global_to_local_interfaces.clear();
    this->local_to_global_interfaces.clear();
}

EndNameSpace
