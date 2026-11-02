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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\\*---------------------------------------------------------------------------*/

#include "Zone.h"
#include "GridHandles.h"
#include <vector>
#include <utility>
#include <memory>
#include "LogFile.h"
#include "ZoneState.h"
#include "Fatal.h"
#include <cstdlib>
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
#include "GridTypes.h"
#include "Task.h"
#include <iostream>
#include <string>


BeginNameSpace( ONEFLOW )


std::vector< Grids > Zone::globalGrids;
int Zone::nLocalZones = 0;
int Zone::flag_test_grid = 0;
std::unique_ptr< InterFaceTopo > Zone::interfaceTopo;

void Zone::ReleaseGrids()
{
    // Destroying each Grids vector releases all owned Grid instances.
    Zone::globalGrids.clear();
    Zone::nLocalZones = 0;
    Zone::flag_test_grid = 0;

    ScalarZone::Reset();
    Zone::interfaceTopo.reset();
}

void Zone::Reset()
{
    Zone::ReleaseGrids();
    ZoneState::Reset();
    GridState::Reset();
}

InterFaceTopo & Zone::GetInterfaceTopo()
{
    if ( ! Zone::interfaceTopo )
    {
        Zone::interfaceTopo = std::make_unique< InterFaceTopo >();
    }
    return *Zone::interfaceTopo;
}

void Zone::InitInterfaceTopo()
{
    Zone::GetInterfaceTopo().InitInterfaceTopo();
}

void Zone::AddGrid( int zid, std::unique_ptr< Grid > grid )
{
    if ( Zone::globalGrids.empty() )
    {
        Zone::globalGrids.resize( static_cast< std::size_t >( ZoneState::nZones ) );
    }
    Zone::globalGrids[ static_cast< std::size_t >( zid ) ].push_back( std::move( grid ) );
}


Grid * Zone::GetGrid( int zid, int gl )
{
    return &Zone::GetGridReference( zid, gl );
}

Grid & Zone::GetGridReference( int zid, int gl )
{
    return GridAt( Zone::globalGrids[ static_cast< std::size_t >( zid ) ], gl );
}

Grid * Zone::GetGrid()
{
    return &Zone::GetGridReference();
}

Grid & Zone::GetGridReference()
{
    return Zone::GetGridReference( ZoneState::zid, GridState::gridLevel );
}

const Grid * Zone::GetConstGrid( int zid, int gl )
{
    return &Zone::GetConstGridReference( zid, gl );
}

const Grid & Zone::GetConstGridReference( int zid, int gl )
{
    return GridAt( Zone::globalGrids[ static_cast< std::size_t >( zid ) ], static_cast< std::size_t >( gl ) );
}

const Grid * Zone::GetConstGrid()
{
    return &Zone::GetConstGridReference();
}

const Grid & Zone::GetConstGridReference()
{
    return Zone::GetConstGridReference( ZoneState::zid, GridState::gridLevel );
}

UnsGrid * Zone::GetUnsGrid()
{
    return ONEFLOW::UnsGridCast( Zone::GetGrid() );
}

const UnsGrid * Zone::GetConstUnsGrid()
{
    return ONEFLOW::UnsGridCast( Zone::GetConstGrid() );
}

Grid * Zone::GetCGrid( Grid * grid )
{
    int level = grid->level + 1;
    int ngrid = static_cast< int >( Zone::globalGrids[ static_cast< std::size_t >( ZoneState::zid ) ].size() );
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
    Zone::InitLayout( fileNameList, std::string() );
}

void Zone::InitLayout(
    StringField & fileNameList,
    const std::string & caseDir )
{
    int nTZones = 0;
    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        std::fstream file;

        if ( caseDir.empty() )
        {
            PIO::OpenPrjFile(
                file, fileNameList[ iFile ],
                std::ios_base::in|std::ios_base::binary );
        }
        else
        {
            PIO::OpenCaseFile(
                file, caseDir, fileNameList[ iFile ],
                std::ios_base::in|std::ios_base::binary );
        }

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

    // Two silent failure modes cost hours of measurement before they were
    // caught, so both are now hard errors with the remedy in the message:
    //   1. a launcher that does not propagate the rank count (srun under PMIx
    //      gave MPI_Comm_size==1 on every rank, so each rank solved all zones
    //      and the run was 10x slower with no diagnostic);
    //   2. a decomposition that cannot be split evenly across the ranks.
    VerifyZoneRankContract();
}

namespace {

int RequestedProcessCountFromEnvironment()
{
    const char * names[] = { "SLURM_NTASKS", "SLURM_NPROCS" };
    for ( const char * name : names )
    {
        const char * value = std::getenv( name );
        if ( value == nullptr || value[ 0 ] == '\0' ) continue;
        const int parsed = std::atoi( value );
        if ( parsed > 0 ) return parsed;
    }
    return 0;
}

std::string CountText( int value )
{
    return std::to_string( value );
}

}

void VerifyZoneRankContract()
{
    const int nZones = ZoneState::nZones;
    const int nProc = Parallel::GetNProc();
    if ( nZones <= 0 || nProc <= 0 ) return;

    // Only a launcher that did not propagate its rank count is a problem.
    // mpirun always sets OMPI_COMM_WORLD_SIZE, so an intentional
    // "mpirun -np 1" inside a larger allocation stays legal.
    const char * launcherSize = std::getenv( "OMPI_COMM_WORLD_SIZE" );
    const bool launcherReportedSize =
        launcherSize != nullptr && launcherSize[ 0 ] != '\0';
    // A Slurm step reports how many tasks it launched.  A direct serial run
    // inside a batch allocation has a single-task step, so the partition step
    // stays legal; an "srun -n N" step reports N and must agree with the MPI
    // size the binary actually got.
    const char * stepTasks = std::getenv( "SLURM_STEP_NUM_TASKS" );
    const bool startedBySlurmStep =
        stepTasks != nullptr && std::atoi( stepTasks ) > 1;
    const int requested = RequestedProcessCountFromEnvironment();
    if ( ! launcherReportedSize && startedBySlurmStep
         && requested > 1 && nProc == 1 )
    {
        Fatal(
            "MPI was initialized with 1 rank while the scheduler requested "
            + CountText( requested ) + ". Every rank would solve all zones "
            "(about " + CountText( requested ) + "x redundant work). Launch "
            "the solver with mpirun -np <ranks> instead of srun, or check "
            "that the MPI launcher matches the MPI library the binary was "
            "built with." );
    }
    if ( nZones % nProc != 0 )
    {
        Fatal(
            "grid has " + CountText( nZones ) + " zones but the run has "
            + CountText( nProc ) + " ranks; the zones cannot be split evenly. "
            "Re-partition the grid for this rank count: run "
            "simutask=Partition with gridObj=Partition, part_uns_file=<single "
            "zone .ofl>, npartproc=" + CountText( nProc ) + " and part_file="
            "<output .ofl>, then solve with gridFileName=<output .ofl>." );
    }
}

void Zone::ReadGrid( StringField & fileNameList )
{
    const std::string sourceCaseDir = GridConfig::GetSourceCaseDir();
    Zone::ReadGrid( fileNameList, sourceCaseDir );
}

void Zone::ReadGrid(
    StringField & fileNameList,
    const std::string & caseDir )
{
    Zone::InitLayout( fileNameList, caseDir );
    int zid = 0;
    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        auto gridGroup = std::make_unique< GridGroup >( zid );

        if ( caseDir.empty() )
        {
            gridGroup->ReadGrid( fileNameList[ iFile ] );
        }
        else
        {
            gridGroup->ReadGrid( fileNameList[ iFile ], caseDir );
        }

        zid += gridGroup->nZones;
    }
    Zone::NormalizeLayout();
}

int Zone::GetNumberOfZoneNeighbors( int zoneId )
{
    return Zone::GetInterfaceTopo().data[ zoneId ].size();
}

int Zone::GetNeighborZoneId( int zoneId, int iNeighbor )
{
    return Zone::GetInterfaceTopo().data[ zoneId ][ iNeighbor ];
}

EndNameSpace
