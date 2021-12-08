/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "SolverDef.h"
#include <iostream>
#include <vector>
#include <algorithm>


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
    int nIFaces = ifaces.size();
    ONEFLOW::HXWrite( databook, this->zonej );
    ONEFLOW::HXWrite( databook, nIFaces );
    ONEFLOW::HXWrite( databook, this->ifaces );
    ONEFLOW::HXWrite( databook, this->recv_ifaces );
}

void ScalarIFaceIJ::ReadInterfaceTopology( DataBook * databook )
{
    ONEFLOW::HXRead( databook, this->zonej );
    int nIFaces = -1;
    ONEFLOW::HXRead( databook, nIFaces );
    this->ifaces.resize( nIFaces );
    this->recv_ifaces.resize( nIFaces );

    ONEFLOW::HXRead( databook, this->ifaces );
    ONEFLOW::HXRead( databook, this->recv_ifaces );
}

ScalarIFace::ScalarIFace()
{
    this->dataSend = new DataStorage();
    this->dataRecv = new DataStorage();
}

ScalarIFace::~ScalarIFace()
{
    delete this->dataSend;
    delete this->dataRecv;
}

void ScalarIFace::AddInterface( int global_interface_id, int neighbor_zoneid, int neighbor_cellid )
{
    int ilocal_interface = this->iglobalfaces.size();
    this->iglobalfaces.push_back( global_interface_id );
    this->zones.push_back( neighbor_zoneid );
    this->cells.push_back( neighbor_cellid );
    this->global_to_local_interfaces[ global_interface_id ] = ilocal_interface;
    this->local_to_global_interfaces[ ilocal_interface ] = global_interface_id;
}

int ScalarIFace::GetLocalInterfaceId( int global_interface_id )
{
    return this->global_to_local_interfaces[ global_interface_id ];
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
    for ( int i = 0; i < globalfaces.size(); ++ i )
    {
        int gid = globalfaces[ i ];
        int lid = this->global_to_local_interfaces[ gid ];
        localfaces.push_back( lid );
    }
    //The neighbor of iZone iNei is jzone, and the jNei neighbor of jZone is iZone
    int jNei = FindINeibor( iZone );
    //std::cout << " zoneid = " << this->zoneid << " iZone() = " << iZone << " jNei = " << jNei << "\n";
    ScalarIFaceIJ & iFaceIJ = this->data[ jNei ];
    iFaceIJ.recv_ifaces = localfaces;
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
    int nSize = zones.size();
    std::set<int> nei_zoneidset;
    for ( int i = 0; i < nSize; ++ i )
    {
        int nei_zoneid = zones[ i ];
        nei_zoneidset.insert( nei_zoneid );
    }

    for ( std::set<int>::iterator iter = nei_zoneidset.begin(); iter != nei_zoneidset.end(); ++ iter )
    {
        ScalarIFaceIJ sij;
        int current_nei_zoneid = * iter;
        //sij.zonei = zoneid;
        sij.zonej = current_nei_zoneid;

        for ( int i = 0; i < nSize; ++ i )
        {
            int nei_zoneid = zones[ i ];
            if ( nei_zoneid == current_nei_zoneid )
            {
                sij.cells.push_back( this->cells[ i ] );
                sij.iglobalfaces.push_back( this->iglobalfaces[ i ] );
                sij.ifaces.push_back( i );
            }
        }
        this->data.push_back( sij );
    }
}

DataStorage * ScalarIFace::GetDataStorage( int iSendRecv )
{
    if ( iSendRecv == SEND_STORAGE )
    {
        return this->dataSend;
    }
    else if ( iSendRecv == RECV_STORAGE )
    {
        return this->dataRecv;
    }
    else
    {
        return 0;
    }
}

void ScalarIFace::WriteInterfaceTopology( DataBook * databook )
{
    int nIFaces = this->GetNIFaces();

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

    std::cout << " nIFaces = " << nIFaces << std::endl;

    if ( nIFaces > 0 )
    {
        this->zones.resize( nIFaces );
        this->target_interfaces.resize( nIFaces );
        this->interface_to_bcface.resize( nIFaces );

        ONEFLOW::HXRead( databook, this->zones               );
        ONEFLOW::HXRead( databook, this->target_interfaces   );
        ONEFLOW::HXRead( databook, this->interface_to_bcface );

        int nNeis = -1;
        ONEFLOW::HXRead( databook, nNeis );
        this->data.resize( nNeis );
        for ( int iNei = 0; iNei < nNeis; ++ iNei )
        {
            ScalarIFaceIJ & iFaceIJ = data[ iNei ];
            iFaceIJ.ReadInterfaceTopology( databook );
        }
    }
}

EndNameSpace
