/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "MultiBlock.h"
#include "Zone.h"
#include "NsCtrl.h"
#include "DataBase.h"
#include "SolverState.h"
#include "SolverDef.h"
#include "SimuDef.h"
#include "WallDist.h"
#include "CmxTask.h"
#include "InterFace.h"
#include "SlipFace.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

MultiBlock::MultiBlock()
{
    ;
}

MultiBlock::~MultiBlock()
{
    ;
}

void MultiBlock::ReadMultiBlockGrid()
{
    StringField gridFileList;

    string gridFileName = ONEFLOW::GetGridFileName();

    gridFileList.push_back( gridFileName );

    Zone::ReadGrid( gridFileList );
}

void MultiBlock::SetUpMultigrid()
{
    SolverState::tid = GRID_SOLVER;
    SsSgTask( "CALC_METRICS" );
}

void MultiBlock::LoadGridAndBuildLink()
{
    MultiBlock::ReadMultiBlockGrid();
    MultiBlock::SetUpMultigrid();
    MultiBlock::InitMultiZoneTopo();
}


void MultiBlock::PrepareFlowGrid()
{
    MultiBlock::LoadGridAndBuildLink();
}

void MultiBlock::ProcessFlowWallDist()
{
    AllocWallDist();

    if ( vis_model.vismodel <= 1 ) return;

    if ( ctrl.startStrategy > 0 )
    {
        LoadWallDist();
    }
    else if ( ctrl.ireadwdst == 0 )
    {
        CreateWallDist();
    }
    else
    {
        LoadWallDist();
    }
}

void MultiBlock::ProcessWallDist()
{
    AllocWallDist();
    CreateWallDist();
}

void CreateWallDist()
{
    SolverState::tid = GRID_SOLVER;
    SsSgTask( "FILL_WALL_STRUCT" );
    SsSgTask( "CALC_WALL_DIST" );
    FreeWallStruct();
    SsSgTask( "WRITE_WALL_DIST" );
}

void LoadWallDist()
{
    SolverState::tid = GRID_SOLVER;
    SsSgTask( "READ_WALL_DIST" );
}

void MultiBlock::AllocWallDist()
{
    SolverState::tid = GRID_SOLVER;
    SsSgTask( "ALLOCATE_WALL_DIST" );
}

void MultiBlock::InitMultiZoneTopo()
{
    ONEFLOW::InitInterfaceTopo();
    ONEFLOW::InitSlipFaceTopo();
    MultiBlock::InitOversetTopo();
}

void MultiBlock::InitOversetTopo()
{
}

string GetGridFileName()
{
    return ONEFLOW::GetDataValue< string >( "gridFileName" );
}

void WalldistSimu()
{
    MultiBlock::LoadGridAndBuildLink();
    MultiBlock::ProcessWallDist();
}


EndNameSpace
