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
// Execute(SimuContext&) receives context; production bodies still call
// the original free functions (ctx unused until a task needs rank/args).

#include "SimuTask.h"
#include "SimuContext.h"
#include "SimuTaskRequire.h"
#include "EulerDomainStateSync.h"
#include "FieldSimu.h"
#include "GridGeneration.h"
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
    void Execute( SimuContext& ctx ) override
    {
        RequireSolveFieldContext( ctx );

        // Single pipeline: stage order lives in FieldPipeline::Run( ctx ).
        FieldSimuRunPipeline( ctx );
    }
};

class CreateGridTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }

    void Execute( SimuContext& ctx ) override
    {
        GenerateGrid( ctx.CaseDir() );
    }
};

class WallDistTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }
    void Execute( SimuContext& /*ctx*/ ) override { WalldistSimu(); }
};

class FunctionTestTask : public ISimuTask
{
public:
    void Execute( SimuContext& /*ctx*/ ) override { FunctionTest(); }
};

class TheoryTask : public ISimuTask
{
public:
    void Execute( SimuContext& /*ctx*/ ) override { TheorySimu(); }
};

class ToyModelTask : public ISimuTask
{
public:
    void Execute( SimuContext& /*ctx*/ ) override { ToyModelSimu(); }
};

class PostTask : public ISimuTask
{
public:
    void Execute( SimuContext& /*ctx*/ ) override { PostSimu(); }
};

class PartitionGridTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }

    // The partition pipeline is already owned by the grid workflow: with
    // gridObj=Partition it reads part_uns_file, METIS-partitions it into
    // npartproc zones and writes part_file.  The task only has to dispatch it.
    void Execute( SimuContext& ctx ) override
    {
        GenerateGrid( ctx.CaseDir() );
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
