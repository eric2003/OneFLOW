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
// Production environment bootstrap for SimuContext.
// Delegates to existing globals (compatibility layer for phase 2).

#include "SimuContext.h"
#include "SimuTask.h"
#include "Prj.h"
#include "ParaFile.h"
#include "Parallel.h"
#include "AccelRuntime.h"
#include "SolverMap.h"
#include "SolverNameList.h"
#include "Zone.h"
#include "ZoneState.h"
#include "GridState.h"
#include "FieldManager.h"
#include "DataBase.h"
#include "HeatFlux.h"
#include "NsCom.h"
#include "NsSolver.h"
#include "INsSolver.h"
#include "TurbSolver.h"
#include "TurbCom.h"
#include "Tolerence.h"
#include "LogFile.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

SimuContext::SimuContext( const std::string& caseDir, bool debug )
    : caseDir_( Prj::ResolveCaseDir( caseDir ) )
{
    Prj::hx_debug = debug;
    Prj::run_from_ide = debug;
}

void SimuContext::ProcessCommandLine()
{
    // Keep the selected case directory as explicit case input.
    const CmdLineOptions opt = Prj::ParseCmdLineArgs( args_ );
    // Store the resolved path in the case context so it remains valid even
    // when a later case updates the legacy Prj static state.
    caseDir_ = Prj::ResolveCaseDir( opt.caseDir );

    // Command-line mode is process-wide; case directory binding belongs to
    // SetupCaseEnvironment() so each case owns its legacy IO binding.
    Prj::hx_debug = opt.debug;
    Prj::run_from_ide = opt.debug;
}

void SimuContext::SetupProcessEnvironment()
{
    std::cout << " OneFLOW is running\n";
    ONEFLOW::SetUpParallelEnvironment();

    rank_ = Parallel::GetPid();
    size_ = Parallel::GetNProc();

    ONEFLOW::InitializeAccelRuntime( rank_, size_ );
    processReady_ = true;
}

void SimuContext::SetupCaseEnvironment()
{
    // Mark case setup as active so exception cleanup can release partial
    // case state if control-file initialization fails halfway through.
    envReady_ = true;

    // Bind legacy case-relative IO to the explicit case before reading it.
    Prj::SetPrjBaseDir( caseDir_ );

    logFile.SetCaseDir( caseDir_ );
    ONEFLOW::ReadControlInfo( caseDir_ );
}

void SimuContext::SetupEnvironment()
{
    SetupProcessEnvironment();
    SetupCaseEnvironment();
}

void SimuContext::TeardownCase()
{
    // Device-backed states must release allocations while the selected
    // accelerator runtime is still alive.
    ClearAccelStates();
    SolverMap::FreeSolverMap();
    SolverNameClass::Reset();
    heat_flux.DeAllocate();
    nscom.Reset();
    NsSolver::Reset();
    INsSolver::Reset();
    TurbSolver::Reset();
    turbcom.Reset();
    Tolerence::Reset();
    Zone::ReleaseGrids();
    ZoneState::Reset();
    GridState::Reset();
    FieldManagerRegistry::FreeFieldManager();
    GetGlobalDataBase()->dataField->Clear();
    GetGlobalDataBase()->dataPara->Clear();
    logFile.ClearCaseDir();
    Prj::ClearPrjBaseDir();
    envReady_ = false;
}

void SimuContext::FinalizeEnvironment()
{
    ONEFLOW::FinalizeAccelRuntime();
    HXFinalize();
    processReady_ = false;
}

void SimuContext::TeardownEnvironment()
{
    TeardownCase();
    FinalizeEnvironment();
}

void SimuContext::ResolveTaskFromControl()
{
    simu_state.Init();
    task_ = simu_state.Task();
    taskName_ = TaskEnumToString( task_ );
    taskResolved_ = true;
}

EndNameSpace
