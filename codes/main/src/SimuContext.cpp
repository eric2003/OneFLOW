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
// Pure / injectable part of SimuContext (no MPI, no control-file I/O).
// Environment bootstrap lives in SimuContextEnv.cpp so unit tests can link
// only this translation unit.

#include "SimuContext.h"
#include "SimuTask.h"
#include "EulerDomainStateLifecycle.h"
#include <stdexcept>
#include <type_traits>

BeginNameSpace( ONEFLOW )

static_assert(
    std::is_nothrow_move_constructible<std::vector<std::string>>::value,
    "SimuContext args must remain nothrow movable");
static_assert(
    std::is_nothrow_move_constructible<std::string>::value,
    "SimuContext strings must remain nothrow movable");
static_assert(
    std::is_nothrow_move_constructible<
        SimuProcessLifecycleDetail::ProcessRuntimeState>::value,
    "Process lifecycle ownership must remain nothrow movable");
static_assert(
    std::is_nothrow_move_constructible<EulerDomainStateRegistry>::value,
    "Accelerator state registry must remain nothrow movable");
static_assert(
    std::is_nothrow_move_constructible<StringField>::value,
    "Expanded solver names must remain nothrow movable");

SimuContext::SimuContext( std::vector<std::string> args )
    : args_( std::move( args ) )
{
}


SimuContext::SimuContext( SimuContext&& other ) noexcept
    : args_( std::move( other.args_ ) )
    , caseDir_( std::move( other.caseDir_ ) )
    , rank_( other.rank_ )
    , size_( other.size_ )
    , task_( other.task_ )
    , taskName_( std::move( other.taskName_ ) )
    , envReady_( other.envReady_ )
    , processLifecycle_( std::move( other.processLifecycle_ ) )
    , taskResolved_( other.taskResolved_ )
    , accelStates_( std::move( other.accelStates_ ) )
    , expandedSolverNames_( std::move( other.expandedSolverNames_ ) )
{
    // A moved-from context must not advertise resources it no longer owns.
    other.args_.clear();
    other.caseDir_.clear();
    other.rank_ = 0;
    other.size_ = 1;
    other.task_ = TaskEnum::SOLVE_FIELD;
    other.taskName_.clear();
    other.envReady_ = false;
    other.taskResolved_ = false;
    other.expandedSolverNames_.clear();
}


void SimuContext::SetParallelInfo( int rank, int size )
{
    if ( size <= 0 )
    {
        throw std::invalid_argument( "SimuContext::SetParallelInfo: size must be > 0" );
    }
    if ( rank < 0 || rank >= size )
    {
        throw std::invalid_argument( "SimuContext::SetParallelInfo: rank out of range" );
    }
    rank_ = rank;
    size_ = size;
}

void SimuContext::SetTask( TaskEnum task )
{
    task_ = task;
    taskName_ = TaskEnumToString( task );
    taskResolved_ = true;
}

void SimuContext::SetTaskByName( const std::string& taskName )
{
    auto it = TaskFilter.find( taskName );
    if ( it == TaskFilter.end() )
    {
        throw std::invalid_argument(
            "SimuContext::SetTaskByName: unknown task \"" + taskName + "\"" );
    }
    task_ = it->second;
    taskName_ = taskName;
    taskResolved_ = true;
}

void SimuContext::MarkEnvironmentReady( bool ready )
{
    envReady_ = ready;
}

EulerDomainState& SimuContext::InitializeAccelState(
    const EulerDomainBackend& backend,
    const EulerDomainProblem& problem,
    const EulerDomainStateKey& key,
    const EulerDomainConstFieldView& field )
{
    return EulerDomainStateLifecycle::Initialize(
        accelStates_, backend, problem, key, field );
}

EulerDomainState& SimuContext::RestartAccelState(
    const EulerDomainBackend& backend,
    const EulerDomainProblem& problem,
    const EulerDomainStateKey& key,
    const EulerDomainConstFieldView& field )
{
    return EulerDomainStateLifecycle::Restart(
        accelStates_, backend, problem, key, field );
}

void SimuContext::SetExpandedSolverNames( const StringField& names )
{
    expandedSolverNames_ = names;
}

void SimuContext::ClearExpandedSolverNames()
{
    expandedSolverNames_.clear();
}

void SimuContext::EnsureExpandedSolverNames( const StringField& names )
{
    if ( !expandedSolverNames_.empty() )
    {
        return;
    }
    expandedSolverNames_ = names;
}


EndNameSpace
