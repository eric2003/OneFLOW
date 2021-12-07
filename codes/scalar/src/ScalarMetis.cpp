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

#include "ScalarMetis.h"
#include "InterFace.h"
#include "CgnsZbase.h"
#include "DataBook.h"
#include "Dimension.h"
#include "FieldPara.h"
#include "Zone.h"
#include "PIO.h"
#include "DataBase.h"
#include "StrUtil.h"
#include "ScalarDataIO.h"
#include "ScalarGrid.h"
#include "MetisGrid.h"
#include "ScalarField.h"
#include "ScalarIFace.h"
#include "ZoneState.h"
#include "ActionState.h"
#include "GridState.h"
#include "ScalarFieldRecord.h"
#include "ScalarAlloc.h"
#include "SolverDef.h"
#include "Prj.h"
#include "FileUtil.h"
#include "Parallel.h"
#include "ScalarZone.h"
#include "HXCgns.h"
#include "HXMath.h"
#include "SmartGrid.h"
#include <iostream>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

ScalarMetis::ScalarMetis()
{
    ;
}

ScalarMetis::~ScalarMetis()
{
    ;
}

void ScalarMetis::Run()
{
    Dim::SetDimension( ONEFLOW::GetDataValue< int >( "dimension" ) );

    vector< ScalarGrid * > input_grids;
    vector< ScalarGrid * > part_grids;

    int dimension = 1;
    std::string root_gridfile = ONEFLOW::GetDataValue< std::string >("root_gridfile");
    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >("scalar_grid_filename");

    int scalar_flag = ONEFLOW::GetDataValue< int >("scalar_flag");

    ScalarReadGrid( root_gridfile, input_grids );
    ScalarGrid * root_grid = input_grids[ 0 ];
    root_grid->CalcMetrics1D();

    int scalar_npart = ONEFLOW::GetDataValue< int >("scalar_npart");
    cout << " scalar_npart = " << scalar_npart << "\n";

    GridPartition gridPartition;
    gridPartition.PartitionGrid( root_grid, scalar_npart, & part_grids );

    ScalarMetisAddZoneGrid( part_grids );
    ScalarDumpGrid( scalar_grid_filename, part_grids );
}

void ScalarMetis::Create1DMesh()
{
    ScalarGrid * grid = new ScalarGrid();

    int scalar_nx = ONEFLOW::GetDataValue< int >("scalar_nx");
    Real scalar_len = ONEFLOW::GetDataValue< int >("scalar_len");

    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >("scalar_grid_filename");

    grid->GenerateGrid( scalar_nx, 0, scalar_len );
    grid->CalcTopology();
    grid->CalcMetrics1D();

    ScalarDumpGrid( scalar_grid_filename, grid );

    delete grid;

    SmartGrid * smart_grid = new SmartGrid();
    smart_grid->Run();
    delete smart_grid;
}

void ScalarMetis::CreateCgnsMesh1D()
{
    ScalarGrid * grid = new ScalarGrid();

    int scalar_nx = ONEFLOW::GetDataValue< int >("scalar_nx");
    Real scalar_len = ONEFLOW::GetDataValue< int >("scalar_len");

    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >("scalar_grid_filename");

    grid->GenerateGrid( scalar_nx, 0, scalar_len );
    grid->CalcTopology();
    grid->CalcMetrics1D();

    ScalarDumpGrid( scalar_grid_filename, grid );

    delete grid;
}

void ScalarMetis::Create1DMeshFromCgns()
{
    ScalarGrid * grid = new ScalarGrid();

    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >("scalar_grid_filename");
    std::string scalar_cgns_filename = ONEFLOW::GetDataValue< std::string >("scalar_cgns_filename");

    std::string cgnsprjFileName = ONEFLOW::GetPrjFileName( scalar_cgns_filename );

    grid->GenerateGridFromCgns( cgnsprjFileName );
    grid->CalcTopology();
    grid->CalcMetrics1D();

    ScalarDumpGrid( scalar_grid_filename, grid );

    delete grid;
}


void ScalarMetisAddZoneGrid( vector< ScalarGrid * > & part_grids )
{
    int nZones = part_grids.size();
    ZoneState::nZones = nZones;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone::AddGrid( iZone, part_grids[ iZone ] );
    }
}

void ScalarReadGrid( const std::string & gridFileName, vector< ScalarGrid * > & grids )
{
    fstream file;
    OpenPrjFile( file, gridFileName, std::ios_base::in|std::ios_base::binary );

    int nZone = -1;

    ONEFLOW::HXRead( & file, nZone );

    ZoneState::pid.resize( nZone );
    ZoneState::zoneType.resize( nZone );

    ONEFLOW::HXRead( & file, ZoneState::pid );
    ONEFLOW::HXRead( & file, ZoneState::zoneType );

    if ( Parallel::zoneMode == 0 )
    {
        for ( int iZone = 0; iZone < nZone; ++ iZone )
        {
            ZoneState::pid[ iZone ] = ( iZone ) % Parallel::nProc;
        }
    }

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        ScalarGrid * grid = new ScalarGrid();
        grid->id = iZone;
        grid->type = ZoneState::zoneType[ iZone ];
        grid->ReadGrid( file );
        grids.push_back( grid );
    }

    ONEFLOW::CloseFile( file );
}

void ScalarDumpGrid( const std::string & gridFileName, ScalarGrid * grid )
{
    vector< ScalarGrid * > grids;
    grids.push_back( grid );
    ScalarDumpGrid( gridFileName, grids );
}

void ScalarDumpGrid( const std::string & gridFileName, vector< ScalarGrid * > & grids )
{
    fstream file;
    OpenPrjFile( file, gridFileName, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc );
    int nZone = static_cast<int>( grids.size() );

    ZoneState::pid.resize( nZone );
    ZoneState::zoneType.resize( nZone );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        ZoneState::pid[ iZone ] = iZone;
        ZoneState::zoneType[ iZone ] = grids[ iZone ]->type;
    }

    ONEFLOW::HXWrite( & file, nZone );
    ONEFLOW::HXWrite( & file, ZoneState::pid );
    ONEFLOW::HXWrite( & file, ZoneState::zoneType );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        grids[ iZone ]->WriteGrid( file );
    }

    ONEFLOW::CloseFile( file );
}


EndNameSpace
