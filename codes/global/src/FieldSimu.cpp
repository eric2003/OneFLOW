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
#include "FieldSimu.h"
#include "SimuContext.h"
#include "EulerDomainStateSync.h"
#include "CmxTaskNames.h"
#include "Iteration.h"
#include "Ctrl.h"
#include "NsCom.h"
#include "UsdData.h"
#include "MultiBlock.h"
#include "SolverMap.h"
#include "SolverCatalog.h"
#include "CmxTask.h"
#include "Multigrid.h"
#include "BcData.h"
#include "GridState.h"
#include <iostream>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

void FieldSimuSetupGlobals()
{
    InitFlowSimuGlobal();
}

void FieldSimuLoadGrid()
{
    MultiBlock::LoadGridAndBuildLink();
}

void FieldSimuPrepareWallDist()
{
    MultiBlock::ProcessFlowWallDist();
}


void FieldSimuCreateSolvers()
{
    // Prefer SolverCatalog name at the pipeline boundary (owns via SolverMap today).
    SolverCatalog::CreateDefault();
}

void FieldSimuCreateSolvers( const SimuContext & ctx )
{
    if ( ctx.HasExpandedSolverNames() )
    {
        SolverCatalog::CreateDefault(
            ONEFLOW::UMESH,
            &ctx.ExpandedSolverNames() );
    }
    else
    {
        SolverCatalog::CreateDefault();
    }
}

void FieldSimuInitFlowField()
{
    // Stage entry: task name enters CmxTask here (no numerical change)
    ONEFLOW::MultiSolverMultiGridTask( kInitFlowFieldTaskName );
}

void FieldSimuRun()
{
    MultigridSolve();
}

void FieldSimuRun( SimuContext & context )
{
    MultigridSolve( context );
}

void FieldPipeline::Run()
{
    FieldSimuSetupGlobals();
    FieldSimuLoadGrid();
    FieldSimuPrepareWallDist();
    FieldSimuCreateSolvers();
    FieldSimuInitFlowField();  // -> MultiSolverMultiGridTask(kInitFlowFieldTaskName)
    FieldSimuRun();
}

void FieldPipeline::Run( SimuContext & ctx )
{
    FieldSimuSetupGlobals();
    FieldSimuLoadGrid();
    FieldSimuPrepareWallDist();
    FieldSimuCreateSolvers( ctx );
    FieldSimuInitFlowField();
    SyncAllEulerDomainStates( ctx );
    FieldSimuRun( ctx );
}

void FieldSimuRunPipeline()
{
    FieldPipeline::Run();
}

void FieldSimuRunPipeline( SimuContext & ctx )
{
    FieldPipeline::Run( ctx );
}

void FieldSimu()
{
    FieldPipeline::Run();
}

void FieldSimu( SimuContext & context )
{
    FieldSimuRunPipeline( context );
}

void InitFlowSimuGlobal()
{
    vis_model.Init();
    ctrl.Init();
    Iteration::Init();
    usd.InitBasic();
}

void InitializeSolver()
{
    // Compatibility alias for older call sites
    FieldSimuInitFlowField();
}

EndNameSpace
