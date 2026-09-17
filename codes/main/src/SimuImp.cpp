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
#include "SimuTask.h"
#include "System.h"
#include "Fatal.h"
#include "SolverNameList.h"
#include "GridState.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

SimuImp::SimuImp( std::vector<std::string>& args )
    : ctx_( std::make_unique<SimuContext>( args ) )
    , args( args )
{
    ctx_->ProcessCommandLine();
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
    ctx_->TeardownEnvironment();
}

void SimuImp::InitSimu()
{
    ctx_->SetupEnvironment();
}

void SimuImp::RunSimu()
{
    ctx_->ResolveTaskFromControl();

    // Production Solve path: same expanded U* names CreateSolvers would have
    // read from SolverNameClass; placed on context so FieldSimuCreateSolvers(ctx)
    // uses the injectable seam. No-op if tests already injected names.
    if ( ctx_->TaskName() == "Solve" )
    {
        ctx_->EnsureExpandedSolverNames(
            SolverNameClass::GetSolverNames( ONEFLOW::UMESH ) );
    }

    auto task = TaskRegistry::Instance().Create( ctx_->TaskName() );
    if ( ! task )
    {
        Fatal( "unknown or unregistered simutask value!!" );
    }

    if ( task->NeedsSystemMap() )
    {
        ConstructSystemMap();
    }

    task->Execute( *ctx_ );
}


EndNameSpace
