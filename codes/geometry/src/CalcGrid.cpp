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


BeginNameSpace( ONEFLOW )

CalcGrid::CalcGrid() = default;

CalcGrid::~CalcGrid() = default;

void CalcGrid::Init( Grids grids )
{
    this->Init( std::move( grids ), GridConfig::FromDataBase() );
}

void CalcGrid::Init( Grids grids, const GridConfig & config )
{
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
    std::fstream file;
    Prj::OpenPrjFile( file, gridFileName, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc );
    const int nZone = GridsSize( grids );

    ZoneState::pid.resize( nZone );
    ZoneState::zoneType.resize( nZone );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        ZoneState::pid[ iZone ] = iZone;
        ZoneState::zoneType[ iZone ] = GridAt( grids, iZone )->type;
    }

    ONEFLOW::HXWrite( & file, nZone );
    ONEFLOW::HXWrite( & file, ZoneState::pid );
    ONEFLOW::HXWrite( & file, ZoneState::zoneType );

    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        std::cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        GridAt( grids, iZone )->WriteGrid( file );
    }

    Prj::CloseFile( file );
}

void CalcGrid::Post()
{
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
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        this->ReconstructLink( iZone );
    }
}

void CalcGrid::ReconstructLink( int iZone )
{
    UnsGrid * grid = static_cast< UnsGrid * >( GridAt( grids, iZone ) );

    InterFace * interFace = grid->interFace.get();
    grid->nIFaces = grid->interFace->nIFaces;

    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nBFaces = grid->nBFaces;
    int nIFaces = interFace->nIFaces;
    int nPBFace = nBFaces - nIFaces;

    IntField & lCell = grid->GetFaceTopo().GetLeftCells();
    IntField & rCell = grid->GetFaceTopo().GetRightCells();

    FacePair facePair;
    for ( int iFace = 0; iFace < nIFaces; ++ iFace )
    {
        int nei_zone_id = interFace->zoneId[ iFace ];
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
            UnsGrid * nei_Grid = static_cast< UnsGrid * >( GridAt( grids, nei_zone_id ) );

            if ( FindMatch( nei_Grid, & facePair ) )
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
    this->iFaceLink->ReconstructInterFace();
}

void CalcGrid::ResetGridScaleAndTranslate()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = GridAt( grids, iZone );
        ONEFLOW::ResetGridScaleAndTranslate( *grid->nodeMesh, this->config );
    }
}

void CalcGrid::GenerateLink()
{
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
        Grid * grid = GridAt( grids, iZone );
        grid->ModifyBcType( BC::NO_BOUNDARY, BC::INTERFACE );
    }
}

void CalcGrid::GenerateLgMapping()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = GridAt( grids, iZone );
        grid->GenerateLgMapping( this->iFaceLink.get() );
    }
}

void CalcGrid::ReGenerateLgMapping()
{
    this->iFaceLink->InitNewLgMapping();

    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = GridAt( grids, iZone );
        grid->ReGenerateLgMapping( this->iFaceLink.get() );
    }

    this->UpdateLgMapping();
    this->UpdateOtherTopologyTerm();
}

void CalcGrid::UpdateLgMapping()
{
    iFaceLink->UpdateLgMapping();
}

void CalcGrid::UpdateOtherTopologyTerm()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = GridAt( grids, iZone );
        grid->UpdateOtherTopologyTerm( this->iFaceLink.get() );
    }
}

void CalcGrid::MatchInterfaceTopology()
{
    const int nZone = GridsSize( grids );
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = GridAt( grids, iZone );
        this->iFaceLink->MatchInterfaceTopology( grid );
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
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        tmp         = yN[ iNode ];
        yN[ iNode ] = zN[ iNode ];
        zN[ iNode ] = - tmp;
    }
}

EndNameSpace
