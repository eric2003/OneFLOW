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

#include "GridGroup.h"
#include "Zone.h"
#include "ZoneState.h"
#include "PIO.h"
#include "Parallel.h"
#include "SolverDef.h"
#include "Prj.h"
#include "BgGrid.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "GridState.h"
#include "ActionState.h"
#include "HXMath.h"
#include "DataBook.h"
#include "Task.h"
#include <iostream>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )


IntField GridGroup::pid;
IntField GridGroup::zoneType;

GridGroup::GridGroup( int zoneStart )
{
    this->zoneStart = zoneStart;
}

GridGroup::~GridGroup()
{
}

void GridGroup::InitZoneLayout( const string & fileName )
{
    fstream file;
    PIO::ParallelOpenPrj( file, fileName, ios_base::in|ios_base::binary );

    this->InitZoneLayout( file );
    this->SetMultiZoneLayout();

    PIO::ParallelClose( file );
}

void GridGroup::InitZoneLayout( fstream & file )
{
    int fid = Parallel::GetFid();

    ONEFLOW::HXReadBcast( file, & nZones, 1, fid );

    pid.resize( nZones );
    zoneType.resize( nZones );

    pid = 0;

    ONEFLOW::HXReadBcast( file, & pid[ 0 ], nZones, fid );

    if ( Parallel::zoneMode == 0 )
    {
        for ( int iZone = 0; iZone < nZones; ++ iZone )
        {
            pid[ iZone ] = ( zoneStart + iZone ) % Parallel::nProc;
        }
    }

    ONEFLOW::HXReadBcast( file, & zoneType[ 0 ], nZones, fid );
}

void GridGroup::SetMultiZoneLayout()
{
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        int ig = zoneStart + iZone;
        ZoneState::pid     [ ig ] = this->pid     [ iZone ];
        ZoneState::zoneType[ ig ] = this->zoneType[ iZone ];
    }
}

void GridGroup::ReadGrid( const string & fileName )
{
    fstream file;

    PIO::ParallelOpenPrj( file, fileName, ios_base::in|ios_base::binary );

    this->InitZoneLayout( file );
    this->SetMultiZoneLayout();

    for ( int iZone = 0; iZone < this->nZones; ++ iZone )
    {
        int zid = this->zoneStart + iZone;
        this->ReadGrid( file, zid );
    }

    PIO::ParallelClose( file );
}

void GridGroup::ReadGrid( fstream & file, int zid )
{
    int spid = 0;
    int rpid = 0;

    Parallel::GetSrPid( zid, spid, rpid );

    if ( Parallel::pid == rpid )
    {
        int gridType = ZoneState::zoneType[ zid ];
        Grid * grid = ONEFLOW::CreateGrid( gridType );
        grid->level = 0;
        grid->id = zid;
        grid->localId = Zone::nLocalZones ++;
        grid->type = gridType;
        Zone::AddGrid( zid, grid );
    }

    DataBook * dataBook = new DataBook();

    ONEFLOW::ReadAbstractData( file, dataBook, spid, rpid );

    ONEFLOW::DataToGrid( dataBook, zid );

    delete dataBook;

}

void ReadAbstractData( fstream & file, DataBook * dataBook, int sendpid, int recvpid, int tag )
{
    if ( Parallel::pid == sendpid )
    {
        dataBook->ReadFile( file );
    }

    dataBook->SendRecv( sendpid, recvpid, tag );
}

void DataToGrid( DataBook * dataBook, int zid )
{
    int spid = 0;
    int rpid = 0;

    Parallel::GetSrPid( zid, spid, rpid );

    if ( Parallel::pid != rpid ) return;

    Grid * grid = Zone::GetGrid( zid, 0 );

    grid->Decode( dataBook );
}

EndNameSpace