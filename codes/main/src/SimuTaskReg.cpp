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
// Behaviour is intentionally identical to the previous switch in SimuImp::RunSimu().

#include "SimuTask.h"
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
    void Execute() override { FieldSimu(); }
};

class CreateGridTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }
    void Execute() override { GenerateGrid(); }
};

class WallDistTask : public ISimuTask
{
public:
    bool NeedsSystemMap() const override { return true; }
    void Execute() override { WalldistSimu(); }
};

class FunctionTestTask : public ISimuTask
{
public:
    void Execute() override { FunctionTest(); }
};

class TheoryTask : public ISimuTask
{
public:
    void Execute() override { TheorySimu(); }
};

class ToyModelTask : public ISimuTask
{
public:
    void Execute() override { ToyModelSimu(); }
};

class PostTask : public ISimuTask
{
public:
    void Execute() override { PostSimu(); }
};

// PARTITION_GRID is present in TaskEnum / TaskFilter but was never handled
// by the old switch. Register a stub that fails clearly if selected.
class PartitionGridTask : public ISimuTask
{
public:
    void Execute() override
    {
        throw std::runtime_error(
            "Task \"Partition\" is registered but not implemented in this build." );
    }
};

// Static registration - same translation unit as the concrete classes.
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
