/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "Zone.h"
#include "LogFile.h"
#include "ZoneState.h"
#include "InterFace.h"
#include "ScalarZone.h"
#include "GridGroup.h"
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


BeginNameSpace( ONEFLOW )


HXVector< Grids * > Zone::globalGrids;
int Zone::nLocalZones = 0;
int Zone::flag_test_grid = 0;

Zone::Zone()
{
}

Zone::~Zone()
{
}

void Zone::AddGrid( int zid, Grid * grid )
{
    if ( Zone::globalGrids.size() == 0 )
    {
        Zone::globalGrids.resize( ZoneState::nZones, 0 );
    }
    Grids * grids = Zone::globalGrids[ zid ];
    if ( ! grids )
    {
        grids = new Grids;
        Zone::globalGrids[ zid ] = grids;
    }
    grids->push_back( grid );
}

Grid * Zone::GetGrid( int zid, int gl )
{
    return ( * Zone::globalGrids[ zid ] )[ gl ];
}

Grid * Zone::GetGrid()
{
    return Zone::GetGrid( ZoneState::zid, GridState::gridLevel );
}

UnsGrid * Zone::GetUnsGrid()
{
    return ONEFLOW::UnsGridCast( Zone::GetGrid() );
}

Grid * Zone::GetCGrid( Grid * grid )
{
    int level = grid->level + 1;
    int ngrid = ( * Zone::globalGrids[ ZoneState::zid ] ).size();
    if ( level >= ngrid ) return 0;
    return Zone::GetGrid( ZoneState::zid, level );
}

Grid * Zone::GetFGrid( Grid * grid )
{
    int level = grid->level - 1;
    level = MAX( level, 0 );
    return Zone::GetGrid( ZoneState::zid, level );
}

void Zone::InitLayout( StringField & fileNameList )
{
    int nTZones = 0;
    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        std::fstream file;
        PIO::OpenPrjFile( file, fileNameList[ iFile ], std::ios_base::in|std::ios_base::binary );

        int nZones = 0;

        ONEFLOW::HXReadBcast( file, & nZones, 1, Parallel::GetFid() );

        nTZones += nZones;

        PIO::CloseFile( file );
    }
    std::cout << " nTZones = " << nTZones << std::endl;
    logFile << "  nTZones = " << nTZones << "\n";

    ZoneState::nZones = nTZones;
    ZoneState::pid.resize( ZoneState::nZones );
    ZoneState::zoneType.resize( ZoneState::nZones );
}

void Zone::NormalizeLayout()
{
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        if ( ! ZoneState::IsValidZone( zId ) ) continue;
        ZoneState::localZid.push_back( zId );
    }
    ZoneState::nLocal = ZoneState::localZid.size();
}

void Zone::ReadGrid( StringField & fileNameList )
{
    Zone::InitLayout( fileNameList );
    int zid = 0;
    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        GridGroup * gridGroup = new GridGroup( zid );
        gridGroup->ReadGrid( fileNameList[ iFile ] );
        zid += gridGroup->nZones;
        delete gridGroup;
    }
    Zone::NormalizeLayout();
}

void Zone::AddScalarGrid( int zid, ScalarGrid * grid )
{
    ScalarZone::AddGrid( zid, grid );
}

ScalarGrid * Zone::GetScalarGrid( int iZone )
{
    return ScalarZone::GetGrid( iZone );
}

ScalarGrid * Zone::GetScalarGrid()
{
    return ScalarZone::GetGrid();
}

int Zone::GetNumberOfZoneNeighbors( int zoneId )
{
    return interFaceTopo.data[ zoneId ].size();
}

int Zone::GetNeighborZoneId( int zoneId, int iNeighbor )
{
    return interFaceTopo.data[ zoneId ][ iNeighbor ];
}

EndNameSpace
