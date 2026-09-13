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
// Concrete ISimuTask wrappers around existing free functions.
// Execute(const SimuContext&) receives context; production bodies still call
// the original free functions (ctx unused until a task needs rank/args).

#include "SimuTask.h"
#include "SimuContext.h"
#include "FieldSimu.h"
#include "GridFactory.h"
#include "MultiBlock.h"
#include "Test.h"
#include "Theory.h"
#include "SimpleSimu.h"
#include "PostProcess.h"
#include <stdexcept>

BeginNameSpace( ONEFLOW )

namespace {

class SolveFieldTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }

    void Execute( const SimuContext& ctx ) override
    {
        // Precondition: environment + control already resolved by SimuImp.
        if ( ! ctx.IsEnvironmentReady() )
        {
            throw std::runtime_error(
                "SolveFieldTask: environment not ready (SetupEnvironment required)" );
        }
        if ( ctx.TaskName() != "Solve" )
        {
            throw std::runtime_error(
                "SolveFieldTask: unexpected task name \"" + ctx.TaskName() + "\"" );
        }

        // Same six stages as FieldSimu(); order must stay identical.
        FieldSimuSetupGlobals();
        FieldSimuLoadGrid();
        FieldSimuPrepareWallDist();
        FieldSimuCreateSolvers();
        FieldSimuInitFlowField();
        FieldSimuRun();
    }
};

class CreateGridTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }
    void Execute( const SimuContext& /*ctx*/ ) override { GenerateGrid(); }
};

class WallDistTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }
    void Execute( const SimuContext& /*ctx*/ ) override { WalldistSimu(); }
};

class FunctionTestTask : public ISimuTask
{
public:
    void Execute( const SimuContext& /*ctx*/ ) override { FunctionTest(); }
};

class TheoryTask : public ISimuTask
{
public:
    void Execute( const SimuContext& /*ctx*/ ) override { TheorySimu(); }
};

class ToyModelTask : public ISimuTask
{
public:
    void Execute( const SimuContext& /*ctx*/ ) override { ToyModelSimu(); }
};

class PostTask : public ISimuTask
{
public:
    void Execute( const SimuContext& /*ctx*/ ) override { PostSimu(); }
};

class PartitionGridTask : public ISimuTask
{
public:
    void Execute( const SimuContext& /*ctx*/ ) override
    {
        throw std::runtime_error(
            "Task \"Partition\" is registered but not implemented in this build." );
    }
};

const bool kTasksRegistered = []() {
    auto& reg = TaskRegistry::Instance();
    reg.Register( "Solve",        []() { return std::make_unique<SolveFieldTask>(); } );
    reg.Register( "Grid",         []() { return std::make_unique<CreateGridTask>(); } );
    reg.Register( "WallDist",     []() { return std::make_unique<WallDistTask>(); } );
    reg.Register( "Partition",    []() { return std::make_unique<PartitionGridTask>(); } );
    reg.Register( "FunctionTest", []() { return std::make_unique<FunctionTestTask>(); } );
    reg.Register( "Theory",       []() { return std::make_unique<TheoryTask>(); } );
    reg.Register( "ToyModel",     []() { return std::make_unique<ToyModelTask>(); } );
    reg.Register( "PostTask",     []() { return std::make_unique<PostTask>(); } );
    return true;
}();

} // anonymous namespace

EndNameSpace
