/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "ZoneState.h"
#include "GridGroup.h"
#include "Parallel.h"
#include "SolverDef.h"
#include "BasicIO.h"
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


HXVector< Grids * > Zone::globalGrids;
int Zone::nLocalZones = 0;

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
    return Zone::GetGrid( ZoneState::zid, level);
}

Grid * Zone::GetFGrid( Grid * grid )
{
    int level = grid->level - 1;
    level = MAX( level, 0 );
    return Zone::GetGrid( ZoneState::zid, level);
}

void Zone::InitLayout( StringField & fileNameList )
{
    int nTZones = 0;
    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        fstream file;
        PIO::ParallelOpenPrj( file, fileNameList[ iFile ], ios_base::in|ios_base::binary );

        int nZones = 0;

        ONEFLOW::HXReadBcast( file, & nZones, 1, Parallel::GetFid() );

        nTZones += nZones;

        PIO::ParallelClose( file );
    }
    cout << " nTZones = " << nTZones << endl;

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

PIO::PIO()
{
    ;
}

PIO::~PIO()
{
    ;
}

string PIO::GetDirName( const string & fileName )
{
    size_t pos = fileName.find_last_of("\\/");
    if ( string::npos == pos )
    {
        return "";
    }
    else
    {
        return fileName.substr(0, pos);
    }
}

void PIO::ParallelOpen( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    PIO::Open( file, fileName, openMode );
}

void PIO::ParallelOpenPrj( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileName;

    string prjFileName = ONEFLOW::StrIO.str();
    string prj_dir = PIO::GetDirName( prjFileName );
    //cout << " prj_dir = " << prj_dir << " prjFileName = " << prjFileName << " fileName = " << fileName << "\n";
    //cout << " PrjStatus::prjBaseDir = " << PrjStatus::prjBaseDir << "\n";
    if ( ! DirExist( prj_dir ) )
    {
        MakeDir( prj_dir );
    }
    PIO::Open( file, prjFileName, openMode );
}

void PIO::ParallelOpenPrj()
{
    PIO::ParallelOpenPrj( * ActionState::file, TaskState::task->fileInfo->fileName, TaskState::task->fileInfo->openMode );
}

void PIO::ParallelClose()
{
    PIO::ParallelClose( * ActionState::file );
}

void PIO::ParallelClose( fstream & file )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    PIO::Close( file );
}

void PIO::Open( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    file.open( fileName.c_str(), openMode );
    if ( ! file )
    {
        cout << "could not open " << fileName << endl;
        Stop("");
    }
}

void PIO::Close( fstream & file )
{
    file.close();
    file.clear();
}


EndNameSpace