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

int ZoneState::nLocal = 0;
int ZoneState::nZones = 0;
int ZoneState::zid = 0;
int ZoneState::szid = 0;
int ZoneState::rzid = 0;

IntField ZoneState::pid;
IntField ZoneState::zoneType;
IntField ZoneState::localZid;

ZoneState::ZoneState()
{
    ;
}

ZoneState::~ZoneState()
{
    ;
}

bool ZoneState::IsValidZone( int zoneId )
{
    return ZoneState::pid[ zoneId ] == Parallel::GetPid();
}

int ZoneState::GetZid( int iSr )
{
    if ( iSr == GREAT_SEND )
    {
        return ZoneState::szid;
    }
    return ZoneState::rzid;
}

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
    cout << " prj_dir = " << prj_dir << " prjFileName = " << prjFileName << " fileName = " << fileName << "\n";
    cout << " PrjStatus::prjBaseDir = " << PrjStatus::prjBaseDir << "\n";
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