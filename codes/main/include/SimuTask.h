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
#pragma once
#include "NamespaceMacros.h"
#include "SimuDef.h"
#include <memory>
#include <functional>
#include <unordered_map>
#include <vector>
#include <string>

BeginNameSpace( ONEFLOW )

class SimuContext; // phase 2+ : tasks receive explicit runtime context

// Strategy / Command interface for a single simulation task.
// Concrete tasks wrap existing free functions so behaviour stays unchanged.
class ISimuTask
{
public:
    virtual ~ISimuTask() = default;

    // Execute with the current run context (rank, args, task name, ...).
    // Production tasks may ignore ctx until they need it; tests can assert on it.
    virtual void Execute( SimuContext& ctx ) = 0;

    // Optional: whether ConstructSystemMap() must run before Execute().
    virtual bool NeedsSystemMap() const { return false; }
};

using SimuTaskCreator = std::function<std::unique_ptr<ISimuTask>()>;

// Registry keyed by the same control-file strings as TaskFilter
// (e.g. "Solve", "Grid", "ToyModel", ...).
class TaskRegistry
{
public:
    static TaskRegistry& Instance();

    void Register( const std::string& taskName, SimuTaskCreator creator );
    std::unique_ptr<ISimuTask> Create( const std::string& taskName ) const;
    std::unique_ptr<ISimuTask> Create( TaskEnum task ) const;
    bool Contains( const std::string& taskName ) const;
    std::vector<std::string> GetAllRegisteredNames() const;

private:
    TaskRegistry() = default;
    std::unordered_map<std::string, SimuTaskCreator> m_creators;
};

// Map TaskEnum back to the canonical control-file string.
const std::string& TaskEnumToString( TaskEnum task );

// Register a concrete task type under a fixed case name.
// Invoke at file scope (after class definition).
#define REGISTER_SIMU_TASK( ClassName, TaskName )                             \
namespace {                                                                   \
bool ClassName##_simu_task_registered = [](){                                 \
    ONEFLOW::TaskRegistry::Instance().Register( TaskName, [](){               \
        return std::make_unique<ClassName>();                                 \
    } );                                                                      \
    return true;                                                              \
}();                                                                          \
}

EndNameSpace
