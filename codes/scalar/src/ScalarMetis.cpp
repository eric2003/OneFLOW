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

#include "ScalarMetis.h"
#include <memory>
#include <utility>
#include "InterFace.h"
#include "CgnsZbase.h"
#include "DataBook.h"
#include "Dimension.h"
#include "FieldPara.h"
#include "Zone.h"
#include "PIO.h"
#include "DataBase.h"
#include "StringUtils.h"
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

#include "Parallel.h"
#include "ScalarZone.h"
#include "HXCgns.h"
#include "HXMath.h"
#include "SmartGrid.h"
#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>


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

    int dimension = 1;
    std::string root_gridfile = ONEFLOW::GetDataValue< std::string >( "root_gridfile" );
    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >( "scalar_grid_filename" );

    int scalar_flag = ONEFLOW::GetDataValue< int >( "scalar_flag" );

    std::vector< std::unique_ptr< ScalarGrid > > input_grids = ScalarReadGrid( root_gridfile );
    ScalarGrid & root_grid = *input_grids[ 0 ];
    root_grid.CalcMetrics1D();

    int scalar_npart = ONEFLOW::GetDataValue< int >( "scalar_npart" );
    std::cout << " scalar_npart = " << scalar_npart << "\n";

    std::vector< std::unique_ptr< ScalarGrid > > part_grids =
        GridPartition::PartitionGrid( root_grid, scalar_npart );

    ScalarDumpGrid( scalar_grid_filename, part_grids );
    ScalarMetisAddZoneGrid( std::move( part_grids ) );
}

void ScalarMetis::Create1DMesh()
{
    auto grid = std::make_unique< ScalarGrid >();

    int scalar_nx = ONEFLOW::GetDataValue< int >( "scalar_nx" );
    Real scalar_len = ONEFLOW::GetDataValue< int >( "scalar_len" );

    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >( "scalar_grid_filename" );

    grid->GenerateGrid( scalar_nx, 0, scalar_len );
    grid->CalcTopology();
    grid->CalcMetrics1D();

    ScalarDumpGrid( scalar_grid_filename, *grid );

    auto smart_grid = std::make_unique< SmartGrid >();
    smart_grid->Run();
}

void ScalarMetis::CreateCgnsMesh1D()
{
    auto grid = std::make_unique< ScalarGrid >();

    int scalar_nx = ONEFLOW::GetDataValue< int >( "scalar_nx" );
    Real scalar_len = ONEFLOW::GetDataValue< int >( "scalar_len" );

    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >( "scalar_grid_filename" );

    grid->GenerateGrid( scalar_nx, 0, scalar_len );
    grid->CalcTopology();
    grid->CalcMetrics1D();

    ScalarDumpGrid( scalar_grid_filename, *grid );
}

void ScalarMetis::Create1DMeshFromCgns()
{
    auto grid = std::make_unique< ScalarGrid >();

    std::string scalar_grid_filename = ONEFLOW::GetDataValue< std::string >( "scalar_grid_filename" );
    std::string scalar_cgns_filename = ONEFLOW::GetDataValue< std::string >( "scalar_cgns_filename" );

    std::string cgnsprjFileName = Prj::GetPrjFileName( scalar_cgns_filename );

    grid->GenerateGridFromCgns( cgnsprjFileName );
    grid->CalcTopology();
    grid->CalcMetrics1D();

    ScalarDumpGrid( scalar_grid_filename, *grid );

}

void ScalarMetisAddZoneGrid( std::vector< std::unique_ptr< ScalarGrid > > part_grids )
{
    int nZones = static_cast< int >( part_grids.size() );
    ZoneState::nZones = nZones;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        ScalarZone::AddGrid( iZone, std::move( part_grids[ iZone ] ) );
    }
}

std::vector< std::unique_ptr< ScalarGrid > > ScalarReadGrid( const std::string & gridFileName )
{
    std::vector< std::unique_ptr< ScalarGrid > > grids;
    std::fstream file;
    Prj::OpenPrjFile( file, gridFileName, std::ios_base::in|std::ios_base::binary );

    int nZone = -1;

    ONEFLOW::HXRead( & file, nZone );
    if ( ! file )
    {
        throw std::runtime_error( "ScalarReadGrid: failed to read the grid zone count" );
    }

    if ( nZone <= 0 )
    {
        Prj::CloseFile( file );
        throw std::runtime_error( "ScalarReadGrid: grid file must contain at least one zone" );
    }

    // Keep file metadata local until every zone has been read successfully.
    // A truncated file must not leave the global runtime layout partially updated.
    IntField zonePids( nZone );
    IntField zoneTypes( nZone );

    ONEFLOW::HXRead( & file, zonePids );
    ONEFLOW::HXRead( & file, zoneTypes );
    if ( ! file )
    {
        throw std::runtime_error( "ScalarReadGrid: truncated zone metadata" );
    }

    if ( Parallel::zoneMode == 0 )
    {
        for ( int iZone = 0; iZone < nZone; ++ iZone )
        {
            zonePids[ iZone ] = iZone % Parallel::nProc;
        }
    }

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        std::cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        auto grid = std::make_unique< ScalarGrid >();
        grid->id = iZone;
        grid->type = zoneTypes[ iZone ];
        grid->ReadGrid( file );
        if ( ! file )
        {
            throw std::runtime_error( "ScalarReadGrid: truncated grid data for zone " + std::to_string( iZone ) );
        }
        grids.push_back( std::move( grid ) );
    }

    // Publish the layout only after the entire file has been parsed.
    ZoneState::pid = std::move( zonePids );
    ZoneState::zoneType = std::move( zoneTypes );

    Prj::CloseFile( file );
    return grids;
}

void ScalarDumpGrid( const std::string & gridFileName, ScalarGrid & grid )
{
    std::fstream file;
    Prj::OpenPrjFile( file, gridFileName, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc );

    const int nZone = 1;
    const std::vector< int > zoneIds = { 0 };
    const std::vector< int > zoneTypes = { grid.type };

    ONEFLOW::HXWrite( & file, nZone );
    ONEFLOW::HXWrite( & file, zoneIds );
    ONEFLOW::HXWrite( & file, zoneTypes );

    std::cout << "iZone = 0 nZone = 1\n";
    grid.WriteGrid( file );

    Prj::CloseOutputFile( file, "ScalarDumpGrid" );
}

void ScalarDumpGrid( const std::string & gridFileName, const std::vector< std::unique_ptr< ScalarGrid > > & grids )
{
    if ( grids.empty() )
    {
        throw std::invalid_argument( "ScalarDumpGrid: at least one grid zone is required" );
    }
    if ( grids.size() > static_cast< size_t >( std::numeric_limits< int >::max() ) )
    {
        throw std::length_error( "ScalarDumpGrid: zone count exceeds the file format limit" );
    }
    for ( const auto & grid : grids )
    {
        if ( ! grid )
        {
            throw std::invalid_argument( "ScalarDumpGrid: grid zone must not be null" );
        }
    }

    // Validate the complete collection before truncating the destination file.
    int nZone = static_cast<int>( grids.size() );
    std::fstream file;
    Prj::OpenPrjFile( file, gridFileName, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc );

    std::vector< int > zoneIds( static_cast< size_t >( nZone ) );
    std::vector< int > zoneTypes( static_cast< size_t >( nZone ) );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        zoneIds[ iZone ] = iZone;
        zoneTypes[ iZone ] = grids[ iZone ]->type;
    }

    ONEFLOW::HXWrite( & file, nZone );
    ONEFLOW::HXWrite( & file, zoneIds );
    ONEFLOW::HXWrite( & file, zoneTypes );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        std::cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        grids[ iZone ]->WriteGrid( file );
    }

    Prj::CloseOutputFile( file, "ScalarDumpGrid" );
}

EndNameSpace
