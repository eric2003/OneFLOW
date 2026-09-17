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
#include "HXDefine.h"
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

// Explicit runtime context for one simulation run (phase 2).
//
// Goal: stop sprinkling Parallel:: / simu_state / accel calls through
// SimuImp.  Internally this type may still *call* the existing globals
// (compatibility layer).  Call sites should depend on SimuContext, not
// on those globals, so they can later be replaced without rewriting
// every task.
class SimuContext
{
public:
    explicit SimuContext( std::vector<std::string> args );
    ~SimuContext() = default;

    SimuContext( const SimuContext& ) = delete;
    SimuContext& operator=( const SimuContext& ) = delete;
    SimuContext( SimuContext&& ) = default;
    SimuContext& operator=( SimuContext&& ) = default;

    // ---- read-only view used by SimuImp / (later) tasks ----
    const std::vector<std::string>& Args() const { return args_; }
    int Rank() const { return rank_; }
    int Size() const { return size_; }
    bool IsEnvironmentReady() const { return envReady_; }
    bool IsTaskResolved() const { return taskResolved_; }
    TaskEnum Task() const { return task_; }
    const std::string& TaskName() const { return taskName_; }

    // Process command line into project globals (existing Prj path).
    void ProcessCommandLine();

    // Production bootstrap: parallel env + control file + accelerator.
    // Still delegates to existing free functions / singletons.
    void SetupEnvironment();

    // Mirror of SetupEnvironment tear-down.
    void TeardownEnvironment();

    // Read simutask from the control database into this context.
    // Requires SetupEnvironment() (or a test SetTask) first for production.
    void ResolveTaskFromControl();

    // ---- test / injection hooks (no MPI, no control file) ----
    void SetParallelInfo( int rank, int size );
    void SetTask( TaskEnum task );
    void SetTaskByName( const std::string& taskName );
    void MarkEnvironmentReady( bool ready = true );

    // Expanded solver registration names (e.g. "UNsSolver").
    // Production Solve path fills these via EnsureExpandedSolverNames before
    // FieldSimuCreateSolvers(ctx); tests may SetExpandedSolverNames directly.
    // Empty means "not yet set" (CreateSolvers falls back to SolverNameClass).
    bool HasExpandedSolverNames() const { return !expandedSolverNames_.empty(); }
    const StringField& ExpandedSolverNames() const { return expandedSolverNames_; }

    // Set expanded names (overwrites). Used by tests and explicit injection.
    void SetExpandedSolverNames( const StringField& names );
    void ClearExpandedSolverNames();

    // Fill only if empty - production preload; does not clobber test injection.
    void EnsureExpandedSolverNames( const StringField& names );
private:
    std::vector<std::string> args_;
    int rank_ = 0;
    int size_ = 1;
    TaskEnum task_ = TaskEnum::SOLVE_FIELD;
    std::string taskName_ = "Solve";
    bool envReady_ = false;
    bool taskResolved_ = false;
    // ...
    StringField expandedSolverNames_;
};

EndNameSpace
