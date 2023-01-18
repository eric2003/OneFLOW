/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

#include "InterFace.h"
#include "ScalarGrid.h"
#include "ScalarIFace.h"
#include "Grid.h"
#include "BgGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "DataStorage.h"
#include "Parallel.h"
#include "SolverDef.h"
#include "LogFile.h"
#include <algorithm>


BeginNameSpace( ONEFLOW )

bool IsValid( InterFace * interFace )
{
    return interFace->nIFaces != 0;
}

InterFace::InterFace()
{
    this->nIFaces = 0;
    this->parent = 0;
    this->AllocSendRecv();
}

InterFace::InterFace( int nIFaces, Grid * parent )
{
    this->AllocSendRecv();
    this->Set( nIFaces, parent );
}

InterFace::~InterFace()
{
    this->DeAllocSendRecv();
    for ( int iNei = 0; iNei < this->interFacePairs.size(); ++ iNei )
    {
        delete this->interFacePairs[ iNei ];
    }
}

void InterFace::AllocSendRecv()
{
    dataSend.resize( MAX_GHOST_LEVELS );
    dataRecv.resize( MAX_GHOST_LEVELS );
    for ( int i = 0; i < MAX_GHOST_LEVELS; ++ i )
    {
        dataSend[ i ] = new DataStorage();
        dataRecv[ i ] = new DataStorage();
    }
}

void InterFace::DeAllocSendRecv()
{
    for ( int i = 0; i < MAX_GHOST_LEVELS; ++ i )
    {
        delete dataSend[ i ];
        delete dataRecv[ i ];
    }
}

void InterFace::Set( int nIFaces, Grid * parent )
{
    this->Resize( nIFaces );
    this->parent = parent;
}

void InterFace::Resize( int nIFaces )
{
    if ( nIFaces <= 0 ) nIFaces = 0;
    this->nIFaces = nIFaces;
    zoneId          .resize( nIFaces );
    localInterfaceId.resize( nIFaces );
    localCellId     .resize( nIFaces );
    i2b .resize( nIFaces );
    idir.resize( nIFaces );
}

void InterFace::InitNeighborFlag( IntField & flags )
{
    flags.resize( ZoneState::nZones, 0 );

    for ( int iFace = 0; iFace < this->nIFaces; ++ iFace )
    {
        flags[ this->zoneId[ iFace ] ] = 1;
    }

    nNeighbor = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        nNeighbor += flags[ iZone ];
    }
    int kkk = 1;
}

void InterFace::AllocateNeighbor()
{
    this->interFacePairs.resize( nNeighbor );

    for ( int iNei = 0; iNei < nNeighbor; ++ iNei )
    {
        this->interFacePairs[ iNei ] = new InterfacePair();
    }
}

void InterFace::InitNeighborZoneInfo()
{
    int nZone = ZoneState::nZones;

    IntField flags;
    this->InitNeighborFlag( flags );

    //The adjacent blocks of this block are calculated
    this->AllocateNeighbor();

    int iNei = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( flags[ iZone ] )
        {
            this->InitNeighborZoneInfo( iNei, iZone );

            ++ iNei;
        }
    }
}

int InterFace::CalcNIFace( int iNei )
{
    int expectedId = this->interFacePairs[ iNei ]->nzid;
    int nIFaceCount = 0;

    for ( int iFace = 0; iFace < this->nIFaces; ++ iFace )
    {
        if ( this->zoneId[ iFace ] == expectedId )
        {
            ++ nIFaceCount;
        }
    }

    return nIFaceCount;
}

void InterFace::InitNeighborZoneInfo( int iNei, int iZone )
{
    InterfacePair * interfacePair = interFacePairs[ iNei ];
    interfacePair->nzid = iZone;

    this->z2n.insert( std::pair< int, int >( iZone, iNei ) );

    int nIFaceCount = this->CalcNIFace( iNei );

    interfacePair->nIFaces = nIFaceCount;
    interfacePair->idsend.resize( nIFaceCount );
    interfacePair->idrecv.resize( nIFaceCount );

    this->FillRecvId( iNei );
}

void InterFace::FillRecvId( int iNei )
{
    InterfacePair * interfacePair = interFacePairs[ iNei ];

    int iCount = 0;
    for ( int iFace = 0; iFace < this->nIFaces; ++ iFace )
    {
        if ( this->zoneId[ iFace ] == interfacePair->nzid )
        {
            //This shows that idrecv is counted locally by the interface of this block
            interfacePair->idrecv[ iCount ] = iFace;
            ++ iCount;
        }
    }
}

void InterFace::CalcSendId( int iNei, IntField & idsend )
{
    InterfacePair * interfacePair = interFacePairs[ iNei ];
    idsend.resize( interfacePair->nIFaces );

    int iCount = 0;
    for ( int iFace = 0; iFace < this->nIFaces; ++ iFace )
    {
        if ( this->zoneId[ iFace ] == interfacePair->nzid )
        {
            //interfacePair->idsend[ iCount ] = this->localInterfaceId[ iFace ];
            idsend[ iCount ] = this->localInterfaceId[ iFace ];
            ++ iCount;
        }
    }
}

void InterFace::SetSendId( int zid, IntField & idsend )
{
    int iNei = this->z2n[ zid ];
    InterfacePair * interfacePair = interFacePairs[ iNei ];
    interfacePair->idsend = idsend;
}

IntField & InterFace::GetInterfaceId( int neiId, int iSr )
{
    if ( iSr == GREAT_SEND )
    {
        return this->interFacePairs[ neiId ]->idsend;
    }
    return this->interFacePairs[ neiId ]->idrecv;
}


InterfacePair::InterfacePair()
{
    ;
}

InterfacePair::~InterfacePair()
{
}

InterFaceTopo::InterFaceTopo()
{
    this->flag_test = 0;
}

InterFaceTopo::~InterFaceTopo()
{
}

void InterFaceTopo::InitInterfaceTopo()
{
    if ( flag_test == 0 )
    {
        this->InitInterfaceTopoImp();
    }
    else
    {
        this->InitInterfaceTopoTest();
    }
}

void InterFaceTopo::InitInterfaceTopoImp()
{
    this->InitZoneNeighborsInfo();
    this->SwapNeighborZoneInfo();
}

void InterFaceTopo::InitInterfaceTopoTest()
{
    this->InitZoneNeighborsInfoTest();
    this->SwapNeighborZoneInfoTest();
}

void InterFaceTopo::InitZoneNeighborsInfoTest()
{
    //This procedure to find all the neighbors of each zone, if the block has docking boundary
    //This neighbor includes the block itself
    //Finally, block numbers of all neighbors are stored in zoneindex in ascending order.
    int nZone = ZoneState::nZones;

    this->data.resize( nZone );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        ScalarGrid * grid = Zone::GetScalarGrid( iZone );

        ScalarIFace * scalarIFace = grid->scalarIFace;

        IntField & neiborZoneIds = this->data[ iZone ];

        int nNei = scalarIFace->data.size();

        for ( int iNei = 0; iNei < nNei; ++ iNei )
        {
            ScalarIFaceIJ & sij = scalarIFace->data[ iNei ];
            neiborZoneIds.push_back( sij.zonej );
        }
    }
}

void InterFaceTopo::SwapNeighborZoneInfoTest()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int pid = ZoneState::pid[ iZone ];

        IntField & neiborZoneIds = this->data[ iZone ];

        int nNeighbor = neiborZoneIds.size();

        ONEFLOW::HXBcast( & nNeighbor, 1, pid );

        if ( nNeighbor == 0 ) continue;

        neiborZoneIds.resize( nNeighbor );
        ONEFLOW::HXBcast( & neiborZoneIds[ 0 ], nNeighbor, pid );
    }
}

void InterFaceTopo::InitZoneNeighborsInfo()
{
    //This procedure to find all the neighbors of each zone, if the block has docking boundary
    //This neighbor includes the block itself
    //Finally, block numbers of all neighbors are stored in zoneindex in ascending order.
    int nZone = ZoneState::nZones;

    this->data.resize( nZone );

    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        Grid * grid = Zone::GetGrid( iZone );

        grid->interFace->InitNeighborZoneInfo();
    }

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        if ( ! ZoneState::IsValidZone( iZone ) ) continue;

        Grid * grid = Zone::GetGrid( iZone );

        IntField & t = this->data[ iZone ];

        for ( int iNei = 0; iNei < grid->interFace->nNeighbor; ++ iNei )
        {
            InterfacePair * interfacePair = grid->interFace->interFacePairs[ iNei ];

            t.push_back( interfacePair->nzid );
        }
        std::sort( t.begin(), t.end() );
    }
}

void InterFaceTopo::SwapNeighborZoneInfo()
{
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int pid = ZoneState::pid[ iZone ];

        IntField & t = this->data[ iZone ];
        int nNeighbor = t.size();

        ONEFLOW::HXBcast( & nNeighbor, 1, pid );

        t.resize( nNeighbor );

        if ( t.size() == 0 ) continue;

        ONEFLOW::HXBcast( & t[ 0 ], nNeighbor, pid );
    }
    this->SwapNeighborsSendContent();
}

void InterFaceTopo::SwapNeighborsSendContent()
{
    int gl = 0;
    for ( int iZone = 0; iZone < ZoneState::nZones; ++ iZone )
    {
        int spid = ZoneState::pid[ iZone ];

        IntField & t = this->data[ iZone ];
        int nNeighbor = t.size();

        for ( int iNei = 0; iNei < nNeighbor; ++ iNei )
        {
            int nZid = t[ iNei ];
            int rpid = ZoneState::pid[ nZid ];
            int nIFaces = 0;
            IntField idsend;

            if ( Parallel::pid == spid )
            {
                Grid * grid = Zone::GetGrid( iZone );
                InterfacePair * interfacePair = grid->interFace->interFacePairs[ iNei ];

                nIFaces = interfacePair->nIFaces;
                
                grid->interFace->CalcSendId( iNei, idsend );
            }

            ONEFLOW::HXSwapData( & nIFaces, 1, spid, rpid, iZone + gl * ZoneState::nZones );

            if ( Parallel::pid  == rpid && 
                          spid  != rpid )
            {
                idsend.resize( nIFaces );
            }

            ONEFLOW::HXSwapData( & idsend[ 0 ], nIFaces, spid, rpid, iZone + gl * ZoneState::nZones );

            if ( Parallel::pid == rpid )
            {
                Grid * gridN = Zone::GetGrid( nZid );
                gridN->interFace->SetSendId( iZone, idsend );
            }
        }
    }
}

InterFaceTopo interFaceTopo;

void InitInterfaceTopo()
{
    interFaceTopo.InitInterfaceTopo();
}

InterFace * InterFaceState::interFace = 0;

InterFaceState::InterFaceState()
{
}

InterFaceState::~InterFaceState()
{
}

EndNameSpace
