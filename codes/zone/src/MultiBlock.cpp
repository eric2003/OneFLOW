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
#include "MultiBlock.h"
#include "Zone.h"
#include "NsCtrl.h"
#include "DataBase.h"
#include "SolverState.h"
#include "SolverDef.h"
#include "SimuDef.h"
#include "WallDist.h"
#include "WallDistPolicy.h"
#include "CmxTask.h"
#include "CmxTaskNames.h"
#include "InterFace.h"
#include "SlipFace.h"
#include <iostream>


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

    // From control database (e.g. cfd.txt: gridFileName = "grid/....ofl")
    std::string gridFileName = ONEFLOW::GetGridFileName();
    gridFileList.push_back( gridFileName );

    // InitLayout (nZones) ¡ú per-file GridGroup::ReadGrid ¡ú NormalizeLayout (localZid)
    Zone::ReadGrid( gridFileList );
}

void MultiBlock::SetUpMultigrid()
{
    SolverState::solverType = GRID_SOLVER;
    SingleSolverSingleGridTask( kCalcMetricsTaskName );
}

void MultiBlock::LoadGridAndBuildLink()
{
    // Stage L1: control DB ¡ú grid path ¡ú Zone layout + binary grid read
    MultiBlock::ReadMultiBlockGrid();

    // Stage L2: geometric metrics (GRID_SOLVER / CALC_METRICS task)
    MultiBlock::SetUpMultigrid();

    // Stage L3: multi-zone topology (interface / slip / overset)
    MultiBlock::InitMultiZoneTopo();
}

void MultiBlock::PrepareFlowGrid()
{
    MultiBlock::LoadGridAndBuildLink();
}

void MultiBlock::ProcessFlowWallDist()
{
    AllocWallDist();

    const FlowWallDistAction action = DecideFlowWallDistAction(
        vis_model.vismodel,
        ctrl.startStrategy,
        ctrl.ireadwdst );

    switch ( action )
    {
    case FlowWallDistAction::SkipAfterAlloc:
        return;
    case FlowWallDistAction::Load:
        LoadWallDist();
        break;
    case FlowWallDistAction::Create:
        CreateWallDist();
        break;
    }
}

void MultiBlock::ProcessWallDist()
{
    AllocWallDist();
    CreateWallDist();
}

void CreateWallDist()
{
    SolverState::solverType = GRID_SOLVER;
    SingleSolverSingleGridTask( kFillWallStructTaskName );
    SingleSolverSingleGridTask( kCalcWallDistTaskName );
    FreeWallStruct();
    SingleSolverSingleGridTask( kWriteWallDistTaskName );
}

void LoadWallDist()
{
    SolverState::solverType = GRID_SOLVER;
    SingleSolverSingleGridTask( kReadWallDistTaskName );
}

void MultiBlock::AllocWallDist()
{
    SolverState::solverType = GRID_SOLVER;
    SingleSolverSingleGridTask( kAllocateWallDistTaskName );
}

void MultiBlock::InitMultiZoneTopo()
{
    ONEFLOW::InitInterfaceTopo();
    ONEFLOW::InitSlipFaceTopo();
    MultiBlock::InitOversetTopo();  // currently empty
}

void MultiBlock::InitOversetTopo()
{
}

std::string GetGridFileName()
{
    return ONEFLOW::GetDataValue< std::string >( "gridFileName" );
}

void WalldistSimu()
{
    MultiBlock::LoadGridAndBuildLink();
    MultiBlock::ProcessWallDist();
}


EndNameSpace
