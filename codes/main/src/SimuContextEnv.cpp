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
// Production environment bootstrap for SimuContext.
// Delegates to existing globals (compatibility layer for phase 2).

#include "SimuContext.h"
#include "SimulationCommandLine.h"
#include "SimuCaseLifecycle.h"
#include "SimuContextParameter.h"
#include "SimuContextIO.h"
#include "SimuTask.h"
#include "Prj.h"
#include "ParaFile.h"
#include "Parallel.h"
#include "AccelRuntime.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include "SolverMap.h"
#include "SolverNameList.h"
#include "Zone.h"
#include "FieldManager.h"
#include "DataBase.h"
#include "DataBaseType.h"
#include "HeatFlux.h"
#include "NsCom.h"
#include "NsSolver.h"
#include "INsSolver.h"
#include "TurbSolver.h"
#include "TurbCom.h"
#include "Tolerance.h"
#include "LogFile.h"
#include "ActionState.h"
#include "Iteration.h"
#include <iostream>
#include <exception>
#include <filesystem>
#include <stdexcept>

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
    const SimulationCommandLineOptions options =
        ParseSimulationCommandLineOptions( args_ );
    // Store the resolved path in the case context so it remains valid even
    // when a later case updates the legacy Prj static state.
    caseDir_ = Prj::ResolveCaseDir( options.caseDirectories.front() );

    // Command-line mode is process-wide; case directory binding belongs to
    // SetupCaseEnvironment() so each case owns its legacy IO binding.
    Prj::hx_debug = options.debug;
    Prj::run_from_ide = options.debug;
}

void SimuContext::SetupProcessEnvironment()
{
    processLifecycle_.Start();
    SimuCaseLifecycleDetail::ExecuteWithRollback(
        [&]()
        {
            std::cout << " OneFLOW is running\n";
            ONEFLOW::SetUpParallelEnvironment();

            rank_ = Parallel::GetPid();
            size_ = Parallel::GetNProc();

            ONEFLOW::InitializeAccelRuntime( rank_, size_ );
        },
        [&]() { FinalizeEnvironment(); } );
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

    // Snapshot and validate all optional roots before applying global path state.
    const SimuContextDetail::CaseIOSettings ioSettings =
        SimuContextDetail::ReadCaseIOSettings(
            RequireGlobalDataBase().RequireDataPara() );
    SimuContextDetail::ApplyCaseIOSettings( ioSettings );
}

void SimuContext::SetupEnvironment()
{
    SetupProcessEnvironment();
    SimuCaseLifecycleDetail::ExecuteWithRollback(
        [&]() { SetupCaseEnvironment(); },
        [&]() { TeardownEnvironment(); } );
}

void SimuContext::TeardownCase()
{
    // Device-backed states and shared HIP solver buffers must release
    // allocations while the selected accelerator runtime is still alive.
    ClearAccelStates();
#ifdef ONEFLOW_ENABLE_HIP
    HipFluxBackend::ReleaseShared();
#endif
    SolverMap::FreeSolverMap();
    SolverNameClass::Reset();
    heat_flux.DeAllocate();
    nscom.Reset();
    NsSolver::Reset();
    INsSolver::Reset();
    TurbSolver::Reset();
    turbcom.Reset();
    Tolerance::Reset();
    Zone::Reset();
    ActionState::Reset();
    Iteration::Reset();
    FieldManagerRegistry::FreeFieldManager();
    DataBase & dataBase = RequireGlobalDataBase();
    dataBase.RequireDataField().Clear();
    dataBase.RequireDataPara().Clear();
    logFile.ClearCaseDir();
    Prj::ClearPrjBaseDir();

    // Clear context-owned case state so the same process can bind a new case
    // without carrying task or solver-name state from the previous case.
    task_ = TaskEnum::SOLVE_FIELD;
    taskName_ = "Solve";
    taskResolved_ = false;
    expandedSolverNames_.clear();
    envReady_ = false;
}

void SimuContext::FinalizeEnvironment()
{
    if ( ! processLifecycle_.IsActive() )
    {
        return;
    }

    std::exception_ptr finalizeFailure;
    try
    {
        ONEFLOW::FinalizeAccelRuntime();
    }
    catch ( ... )
    {
        finalizeFailure = std::current_exception();
    }

    try
    {
        HXFinalize();
    }
    catch ( ... )
    {
        if ( ! finalizeFailure )
        {
            finalizeFailure = std::current_exception();
        }
    }

    processLifecycle_.MarkFinalized();
    if ( finalizeFailure )
    {
        std::rethrow_exception( finalizeFailure );
    }
}

void SimuContext::TeardownEnvironment()
{
    SimuCaseLifecycleDetail::ExecuteWithCleanup(
        [&]() { TeardownCase(); },
        [&]() { FinalizeEnvironment(); } );
}

void SimuContext::ResolveTaskFromControl()
{
    simu_state.Init();
    task_ = simu_state.Task();
    taskName_ = TaskEnumToString( task_ );
    taskResolved_ = true;
}

EndNameSpace
