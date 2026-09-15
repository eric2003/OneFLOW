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
#include <stdexcept>

BeginNameSpace( ONEFLOW )

SimuContext::SimuContext( std::vector<std::string> args )
    : args_( std::move( args ) )
{
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

void SimuContext::SetExpandedSolverNames( const StringField& names )
{
    expandedSolverNames_ = names;
}

void SimuContext::ClearExpandedSolverNames()
{
    expandedSolverNames_.clear();
}

EndNameSpace
