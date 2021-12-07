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

#include "GridGroup.h"
#include "Zone.h"
#include "ZoneState.h"
#include "ScalarGrid.h"

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


GridGroup::GridGroup( int zoneStart )
{
    this->zoneStart = zoneStart;
}

GridGroup::~GridGroup()
{
}

void GridGroup::InitZoneLayout( const std::string & fileName )
{
    std::fstream file;
    PIO::OpenPrjFile( file, fileName, std::ios_base::in|std::ios_base::binary );

    this->InitZoneLayout( file );
    this->SetMultiZoneLayout();

    PIO::CloseFile( file );
}

void GridGroup::InitZoneLayout( std::fstream & file )
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

void GridGroup::ReadGrid( const std::string & fileName )
{
    std::fstream file;

    PIO::OpenPrjFile( file, fileName, std::ios_base::in|std::ios_base::binary );

    this->InitZoneLayout( file );
    this->SetMultiZoneLayout();

    for ( int iZone = 0; iZone < this->nZones; ++ iZone )
    {
        int zid = this->zoneStart + iZone;
        this->ReadGrid( file, zid );
    }

    PIO::CloseFile( file );
}

void GridGroup::ReadGrid( std::fstream & file, int zid )
{
    int spid = 0;
    int rpid = 0;

    Parallel::GetSrPid( zid, spid, rpid );

    if ( Parallel::pid == rpid )
    {
        this->CreateGrid( zid );
    }

    DataBook * dataBook = new DataBook();

    ONEFLOW::ReadAbstractData( file, dataBook, spid, rpid );

    ONEFLOW::DataToGrid( dataBook, zid );

    delete dataBook;
}

void GridGroup::CreateGrid( int zoneId )
{
    if ( Zone::flag_test_grid == 0 )
    {
        this->CreateGridImp( zoneId );
    }
    else
    {
        this->CreateGridTest( zoneId );
    }
}

void GridGroup::CreateGridImp( int zoneId )
{
    int gridType = ZoneState::zoneType[ zoneId ];
    Grid * grid = ONEFLOW::CreateGrid( gridType );
    grid->level = 0;
    grid->id = zoneId;
    grid->localId = Zone::nLocalZones ++;
    grid->type = gridType;
    Zone::AddGrid( zoneId, grid );
}

void GridGroup::CreateGridTest( int zoneId )
{
    int gridType = ZoneState::zoneType[ zoneId ];

    ScalarGrid * grid = new ScalarGrid();
    grid->level = 0;
    grid->id = zoneId;
    grid->localId = Zone::nLocalZones ++;
    grid->type = gridType;

    Zone::AddScalarGrid( zoneId, grid );
}

void ReadAbstractData( std::fstream & file, DataBook * dataBook, int sendpid, int recvpid, int tag )
{
    if ( Parallel::pid == sendpid )
    {
        dataBook->ReadFile( file );
    }

    dataBook->SendRecv( sendpid, recvpid, tag );
}

void DataToGrid( DataBook * dataBook, int zid )
{
    if ( Zone::flag_test_grid == 0 )
    {
        DataToGridImp( dataBook, zid );
    }
    else
    {
        DataToGridTest( dataBook, zid );
    }
}

void DataToGridImp( DataBook * dataBook, int zid )
{
    int spid = 0;
    int rpid = 0;

    Parallel::GetSrPid( zid, spid, rpid );

    if ( Parallel::pid != rpid ) return;

    Grid * grid = Zone::GetGrid( zid, 0 );

    grid->Decode( dataBook );
}

void DataToGridTest( DataBook * dataBook, int zid )
{
    int spid = 0;
    int rpid = 0;

    Parallel::GetSrPid( zid, spid, rpid );

    if ( Parallel::pid != rpid ) return;

    ScalarGrid * grid = Zone::GetScalarGrid( zid );
    grid->ReadGrid( dataBook );
}

EndNameSpace
