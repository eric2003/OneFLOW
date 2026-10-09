/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) under the terms of the GNU General Public License
    as published by the Free Software Foundation.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "CalcGrid.h"
#include "UnsGrid.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "GridTypes.h"
#include "LogFile.h"
#include "HXMath.h"
#include "IFaceLink.h"
#include "InterFace.h"
#include "FaceTopo.h"
#include "Zone.h"
#include "ZoneState.h"
#include "Partition.h"
#include "DataBase.h"
#include "DataBaseIO.h"
#include "Boundary.h"

#include "Fatal.h"
#include "Prj.h"
#include <iostream>
#include <utility>
#include <stdexcept>
#include <limits>


BeginNameSpace( ONEFLOW )

namespace
{
void ValidateGridCollection( const Grids & grids )
{
    if ( grids.empty() )
    {
        throw std::invalid_argument( "CalcGrid: at least one grid zone is required" );
    }
    if ( grids.size() > static_cast< size_t >( std::numeric_limits< int >::max() ) )
    {
        throw std::length_error( "CalcGrid: zone count exceeds the file format limit" );
    }
    for ( const auto & grid : grids )
    {
        if ( ! grid )
        {
            throw std::invalid_argument( "CalcGrid: grid zone must not be null" );
        }
    }
}

UnsGrid & RequireUnsGrid( Grid & grid, const char * errorMessage )
{
    UnsGrid * unsGrid = dynamic_cast< UnsGrid * >( &grid );
    if ( ! unsGrid )
    {
        throw std::invalid_argument( errorMessage );
    }
    return *unsGrid;
}
} // namespace

CalcGrid::CalcGrid() = default;

CalcGrid::~CalcGrid() = default;

IFaceLink & CalcGrid::GetInterfaceLink()
{
    if ( ! this->iFaceLink )
    {
        throw std::logic_error( "CalcGrid: interface link has not been generated" );
    }
    return *this->iFaceLink;
}

void CalcGrid::Init( Grids grids )
{
    this->Init( std::move( grids ), GridConfig::FromDataBase() );
}

void CalcGrid::Init( Grids grids, const GridConfig & config )
{
    // Validate before taking ownership; Post() traverses this collection.
    ValidateGridCollection( grids );

    this->grids = std::move( grids );
    this->config = config;

    if ( this->config.objective == GridObjective::Partition )
    {
        this->gridFileName = this->config.partitionFile;
    }
    else
    {
        this->gridFileName = this->config.targetFile;
    }
}

void CalcGrid::BuildInterfaceLink()
{
    // This public operation can be called independently of Post().
    ValidateGridCollection( grids );

    if ( this->config.objective == GridObjective::Partition )
    {
        const int partitionType = this->config.partitionType;
        if ( partitionType == 1 )
        {
            this->ReconstructLink();
        }
        else
        {
            this->GenerateLink();
        }
    }
    else
    {
        this->GenerateLink();
    }
}

void CalcGrid::Dump()
{
    // Public collection state can be changed after Init(), so validate again
    // before traversing it or truncating the destination file.
    ValidateGridCollection( grids );

    // Validate the complete collection before truncating the destination file.
    const int nZone = static_cast< int >( grids.size() );
    IntField zonePids( nZone );
    IntField zoneTypes( nZone );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        zonePids[ iZone ] = iZone;
        zoneTypes[ iZone ] = GridAt( grids, iZone ).type;
    }

    std::fstream file;
    Prj::OpenPrjFile( file, gridFileName, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc );

    ONEFLOW::HXWrite( & file, nZone );
    ONEFLOW::HXWrite( & file, zonePids );
    ONEFLOW::HXWrite( & file, zoneTypes );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        std::cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        GridAt( grids, iZone ).WriteGrid( file );
    }

    file.flush();
    if ( ! file )
    {
        throw std::runtime_error( "CalcGrid::Dump: failed to write grid data" );
    }
    Prj::CloseFile( file );

    // Publish zone metadata only after the grid output succeeds.
    ZoneState::pid = std::move( zonePids );
    ZoneState::zoneType = std::move( zoneTypes );
}

void CalcGrid::Post()
{
    // Validate the mutable public collection before any post-processing side effects.
    ValidateGridCollection( grids );

    logFile << "GenerateOverset\n";
    this->GenerateOverset();
    logFile << "BuildInterfaceLink\n";
    this->BuildInterfaceLink();
    logFile << "ResetGridScaleAndTranslate\n";
    this->ResetGridScaleAndTranslate();
    logFile << "CalcGrid::Post() Final \n";
}

void CalcGrid::GenerateOverset()
{
}

void CalcGrid::ReconstructLink()
{
    // This public operation may be called without BuildInterfaceLink().
    ValidateGridCollection( grids );

    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        this->ReconstructLink( iZone );
    }
}

void CalcGrid::ReconstructLink( int iZone )
{
    if ( iZone < 0 || static_cast< size_t >( iZone ) >= grids.size() )
    {
        throw std::out_of_range( "CalcGrid::ReconstructLink: zone index is out of range" );
    }

    Grid & baseGrid = GridAt( grids, static_cast< size_t >( iZone ) );
    UnsGrid & grid = RequireUnsGrid(
        baseGrid, "CalcGrid::ReconstructLink: zone must use an unstructured grid" );

    InterFace * interFace = grid.interFace.get();

    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    grid.nIFaces = interFace->nIFaces;

    int nBFaces = grid.nBFaces;
    int nIFaces = interFace->nIFaces;
    int nPBFace = nBFaces - nIFaces;

    IntField & lCell = grid.GetFaceTopo().GetLeftCells();
    IntField & rCell = grid.GetFaceTopo().GetRightCells();

    FacePair facePair;
    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
    {
        int nei_zone_id = interFace->zoneId[ iFace ];
        if ( nei_zone_id < 0 || static_cast< size_t >( nei_zone_id ) >= grids.size() )
        {
            throw std::out_of_range( "CalcGrid::ReconstructLink: interface neighbor zone is out of range" );
        }

        int lc = lCell[ iFace + nPBFace ];
        int rc = rCell[ iFace + nPBFace ];
        int cellIndex  = MAX( lc, rc );
        facePair.lf.zone_id = iZone;
        facePair.lf.face_id = iFace;
        facePair.lf.cell_id = cellIndex;

        facePair.rf.zone_id = nei_zone_id;
        facePair.rf.cell_id = interFace->localCellId[ iFace ];

        if ( nei_zone_id >= iZone )
        {
            Grid & neighborBaseGrid = GridAt( grids, static_cast< size_t >( nei_zone_id ) );
            UnsGrid & neiGrid = RequireUnsGrid(
                neighborBaseGrid,
                "CalcGrid::ReconstructLink: interface neighbor zone must use an unstructured grid" );

            if ( FindMatch( neiGrid, facePair ) )
            {
                interFace->localInterfaceId[ iFace ] = facePair.rf.face_id;
            }
            else
            {
                Fatal("");
            }
        }
    }
}

void CalcGrid::ReconstructInterFace()
{
    this->GetInterfaceLink().ReconstructInterFace();
}

void CalcGrid::ResetGridScaleAndTranslate()
{
    // Check all node meshes before transforming any zone, so a malformed
    // collection cannot leave earlier zones scaled while later zones fail.
    ValidateGridCollection( grids );
    for ( const auto & grid : grids )
    {
        if ( ! grid->nodeMesh )
        {
            throw std::logic_error(
                "CalcGrid::ResetGridScaleAndTranslate: grid node mesh is not initialized" );
        }
    }

    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid & grid = GridAt( grids, iZone );
        ONEFLOW::ResetGridScaleAndTranslate( *grid.nodeMesh, this->config );
    }
}

void CalcGrid::GenerateLink()
{
    // Validate before constructing an interface-link object from the collection.
    ValidateGridCollection( grids );

    this->iFaceLink = std::make_unique< IFaceLink >( grids );

    this->ModifyBcType();

    this->GenerateLgMapping();

    this->ReconstructInterFace();

    this->ReGenerateLgMapping();

    this->MatchInterfaceTopology();
}

void CalcGrid::ModifyBcType()
{
    if ( this->config.ignoreNoBoundary ) return;

    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid & grid = GridAt( grids, iZone );
        grid.ModifyBcType( BC::NO_BOUNDARY, BC::INTERFACE );
    }
}

void CalcGrid::GenerateLgMapping()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid & grid = GridAt( grids, iZone );
        grid.GenerateLgMapping( this->GetInterfaceLink() );
    }
}

void CalcGrid::ReGenerateLgMapping()
{
    this->GetInterfaceLink().InitNewLgMapping();

    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid & grid = GridAt( grids, iZone );
        grid.ReGenerateLgMapping( this->GetInterfaceLink() );
    }

    this->UpdateLgMapping();
    this->UpdateOtherTopologyTerm();
}

void CalcGrid::UpdateLgMapping()
{
    this->GetInterfaceLink().UpdateLgMapping();
}

void CalcGrid::UpdateOtherTopologyTerm()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid & grid = GridAt( grids, iZone );
        grid.UpdateOtherTopologyTerm( this->GetInterfaceLink() );
    }
}

void CalcGrid::MatchInterfaceTopology()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid & grid = GridAt( grids, iZone );
        this->GetInterfaceLink().MatchInterfaceTopology( grid );
    }
}

void CalcGrid::GenerateMultiZoneCalcGrids( Grids grids )
{
    this->GenerateMultiZoneCalcGrids(
        std::move( grids ), GridConfig::FromDataBase() );
}

void CalcGrid::GenerateMultiZoneCalcGrids(
    Grids grids,
    const GridConfig & config )
{
    RegionNameMap::DumpRegion();

    this->Init( std::move( grids ), config );
    this->Post();
    this->Dump();
}

int GetIgnoreNoBc()
{
    return GridConfig::FromDataBase().ignoreNoBoundary ? 1 : 0;
}

std::string GetTargetGridFileName()
{
    return GridConfig::FromDataBase().targetFile;
}

void GenerateMultiZoneCalcGrids(
    Grids grids,
    const GridConfig & config )
{
    CalcGrid calcGrid;
    calcGrid.GenerateMultiZoneCalcGrids( std::move( grids ), config );
}

void GenerateMultiZoneCalcGrids( Grids grids )
{
    GenerateMultiZoneCalcGrids(
        std::move( grids ), GridConfig::FromDataBase() );
}


void ResetGridScaleAndTranslate( NodeMesh & nodeMesh, const GridConfig & config )
{
    const Real scale = config.scale;
    const auto & translate = config.translate;

    const size_t nNodes = nodeMesh.GetNumberOfNodes();

    for ( size_t iNode = 0; iNode < nNodes; ++ iNode )
    {
        nodeMesh.xN[ iNode ] *= scale;
        nodeMesh.yN[ iNode ] *= scale;
        nodeMesh.zN[ iNode ] *= scale;

        nodeMesh.xN[ iNode ] += translate[ 0 ];
        nodeMesh.yN[ iNode ] += translate[ 1 ];
        nodeMesh.zN[ iNode ] += translate[ 2 ];
    }

    if ( config.axisDirection == GridAxisDirection::ZToY )
    {
        TurnZAxisToYAxis( nodeMesh );
    }
}

void ResetGridScaleAndTranslate( NodeMesh & nodeMesh )
{
    ResetGridScaleAndTranslate( nodeMesh, GridConfig::FromDataBase() );
}

void TurnZAxisToYAxis( NodeMesh & nodeMesh )
{
    size_t nNodes = nodeMesh.GetNumberOfNodes();

    RealField & xN = nodeMesh.xN;
    RealField & yN = nodeMesh.yN;
    RealField & zN = nodeMesh.zN;

    Real tmp;
    for ( size_t iNode = 0; iNode < nNodes; ++ iNode )
    {
        tmp         = yN[ iNode ];
        yN[ iNode ] = zN[ iNode ];
        zN[ iNode ] = - tmp;
    }
}

EndNameSpace
