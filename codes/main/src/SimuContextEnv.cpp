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
#include <iostream>

BeginNameSpace( ONEFLOW )

void SimuContext::ProcessCommandLine()
{
    Prj::ProcessCmdLineArgs( args_ );
}

void SimuContext::SetupEnvironment()
{
    std::cout << " OneFLOW is running\n";
    ONEFLOW::SetUpParallelEnvironment();
    ONEFLOW::ReadControlInfo();

    rank_ = Parallel::GetPid();
    size_ = Parallel::GetNProc();

    ONEFLOW::InitializeAccelRuntime( rank_, size_ );
    envReady_ = true;
}

void SimuContext::TeardownEnvironment()
{
    // Device-backed states must release allocations while the selected
    // accelerator runtime is still alive.
    ClearAccelStates();
    ONEFLOW::FinalizeAccelRuntime();
    HXFinalize();
    envReady_ = false;
}

void SimuContext::ResolveTaskFromControl()
{
    simu_state.Init();
    task_ = simu_state.Task();
    taskName_ = TaskEnumToString( task_ );
    taskResolved_ = true;
}

EndNameSpace
