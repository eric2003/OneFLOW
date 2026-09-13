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
#include "SimuImp.h"
#include "SimuDef.h"
#include "SimuTask.h"
#include "System.h"
#include "Prj.h"
#include "ParaFile.h"
#include "Parallel.h"
#include "Fatal.h"
#include "AccelRuntime.h"
#include <iostream>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

SimuImp::SimuImp( std::vector<std::string>& args )
{
    this->args = args;
    Prj::ProcessCmdLineArgs( args );
}

SimuImp::~SimuImp()
{
}

void SimuImp::Run()
{
    this->PreProcess();
    this->MainProcess();
    this->PostProcess();
}

void SimuImp::PreProcess()
{
    InitSimu();
}

void SimuImp::MainProcess()
{
    RunSimu();
}

void SimuImp::PostProcess()
{
    ONEFLOW::FinalizeAccelRuntime();
    HXFinalize();
}

void SimuImp::RunSimu()
{
    // Resolve task type from control database (unchanged behaviour).
    simu_state.Init();

    const TaskEnum taskEnum = simu_state.Task();
    const std::string& taskName = TaskEnumToString( taskEnum );

    auto task = TaskRegistry::Instance().Create( taskName );
    if ( ! task )
    {
        Fatal( "unknown or unregistered simutask value!!" );
    }

    // Preserve the previous ConstructSystemMap gate.
    if ( task->NeedsSystemMap() )
    {
        ConstructSystemMap();
    }

    task->Execute();
}

void SimuImp::InitSimu()
{
    std::cout << " OneFLOW is running\n";
    ONEFLOW::SetUpParallelEnvironment();
    ONEFLOW::ReadControlInfo();
    ONEFLOW::InitializeAccelRuntime(
        Parallel::GetPid(),
        Parallel::GetNProc() );
}


EndNameSpace
