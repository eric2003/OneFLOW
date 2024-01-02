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

#include "CalcGrid.h"
#include "UnsGrid.h"
#include "Grid.h"
#include "GridPara.h"
#include "NodeMesh.h"
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

#include "Stop.h"
#include "Prj.h"
#include "HXPointer.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

CalcGrid::CalcGrid()
{
    iFaceLink = 0;
}

CalcGrid::~CalcGrid()
{
    delete iFaceLink;
}

void CalcGrid::Init( Grids & grids )
{
    this->grids = grids;
    this->grids.SetDeleteFlag( true );
    int gridObj = GetDataValue< int >( "gridObj" );
    if ( gridObj == 3 )
    {
        std::string part_uns_file = GetDataValue< std::string >( "part_uns_file" );
        this->gridFileName = part_uns_file;
    }
    else
    {
        this->gridFileName = ONEFLOW::GetTargetGridFileName();
    }
}

void CalcGrid::BuildInterfaceLink()
{
    int gridObj = GetDataValue< int >( "gridObj" );

    if ( gridObj == 3 )
    {
        int partition_type = GetDataValue< int >( "partition_type" );
        if ( partition_type == 1 )
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
    //std::cout << __FUNCTION__ << std::endl;
    std::fstream file;
    Prj::OpenPrjFile( file, gridFileName, std::ios_base::out|std::ios_base::binary|std::ios_base::trunc );
    int nZone = static_cast<int>(grids.size());

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
        std::cout << "iZone = " << iZone << " nZone = " << nZone << "\n";
        grids[ iZone ]->WriteGrid( file );
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
    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        this->ReconstructLink( iZone );
    }
}

void CalcGrid::ReconstructLink( int iZone )
{
    UnsGrid * grid = UnsGridCast( grids[ iZone ] );

    InterFace * interFace = grid->interFace;
    grid->nIFaces = grid->interFace->nIFaces;

    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nBFaces = grid->nBFaces;
    int nIFaces = interFace->nIFaces;
    int nPBFace = nBFaces - nIFaces;

    IntField & lCell = grid->faceTopo->lCells;
    IntField & rCell = grid->faceTopo->rCells;

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
            UnsGrid * nei_Grid = UnsGridCast( grids[ nei_zone_id ] );

            if ( FindMatch( nei_Grid, & facePair ) )
            {
                interFace->localInterfaceId[ iFace ] = facePair.rf.face_id;
            }
            else
            {
                Stop("");
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
    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = grids[ iZone ];
        ONEFLOW::ResetGridScaleAndTranslate( grid->nodeMesh );
    }
}

void CalcGrid::GenerateLink()
{
    this->iFaceLink = new IFaceLink( grids );

    this->ModifyBcType();

    this->GenerateLgMapping();

    this->ReconstructInterFace();

    this->ReGenerateLgMapping();

    this->MatchInterfaceTopology();
}

void CalcGrid::ModifyBcType()
{
    int ignoreNoBc = ONEFLOW::GetIgnoreNoBc();

    if ( ignoreNoBc ) return;

    //change NO_BOUNDARY to INTERFACE
    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = grids[ iZone ];
        grid->ModifyBcType( BC::NO_BOUNDARY, BC::INTERFACE );
    }
}

void CalcGrid::GenerateLgMapping()
{
    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = grids[ iZone ];
        grid->GenerateLgMapping( this->iFaceLink );
    }
}

void CalcGrid::ReGenerateLgMapping()
{
    this->iFaceLink->InitNewLgMapping();

    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = grids[ iZone ];
        grid->ReGenerateLgMapping( this->iFaceLink );
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
    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = grids[ iZone ];
        grid->UpdateOtherTopologyTerm( this->iFaceLink );
    }
}

void CalcGrid::MatchInterfaceTopology()
{
    int nZone = static_cast<int>(grids.size());
    for ( int iZone = 0; iZone < nZone; ++ iZone )
    {
        Grid * grid = grids[ iZone ];
        this->iFaceLink->MatchInterfaceTopology( grid );
    }
}

void CalcGrid::GenerateMultiZoneCalcGrids( Grids & grids )
{
    RegionNameMap::DumpRegion();

    this->Init( grids );
    this->Post();
    this->Dump();
}

int GetIgnoreNoBc()
{
    return ONEFLOW::GetDataValue< int >( "ignoreNoBc" );
}

std::string GetTargetGridFileName()
{
    return ONEFLOW::GetDataValue< std::string >( "targetGridFileName" );
}

void GenerateMultiZoneCalcGrids( Grids & grids )
{
    CalcGrid * calcGrid = new CalcGrid();
    calcGrid->GenerateMultiZoneCalcGrids( grids );
    delete calcGrid;
}

void ResetGridScaleAndTranslate( NodeMesh * nodeMesh )
{
    size_t nNodes = nodeMesh->GetNumberOfNodes();

    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        nodeMesh->xN[ iNode ] *= grid_para.gridScale;
        nodeMesh->yN[ iNode ] *= grid_para.gridScale;
        nodeMesh->zN[ iNode ] *= grid_para.gridScale;

        nodeMesh->xN[ iNode ] += grid_para.gridTrans[ 0 ];
        nodeMesh->yN[ iNode ] += grid_para.gridTrans[ 1 ];
        nodeMesh->zN[ iNode ] += grid_para.gridTrans[ 2 ];
    }

    if ( grid_para.axis_dir == 1 )
    {
        TurnZAxisToYAxis( nodeMesh );
    }
}

void TurnZAxisToYAxis( NodeMesh * nodeMesh )
{
    size_t nNodes = nodeMesh->GetNumberOfNodes();

    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    Real tmp;
    for ( int iNode = 0; iNode < nNodes; ++ iNode )
    {
        tmp         = yN[ iNode ];
        yN[ iNode ] = zN[ iNode ];
        zN[ iNode ] = - tmp;
    }
}

EndNameSpace
