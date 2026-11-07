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
#include "SimuTask.h"

BeginNameSpace( ONEFLOW )

TaskRegistry& TaskRegistry::Instance()
{
    static TaskRegistry inst;
    return inst;
}

void TaskRegistry::Register( const std::string& taskName, SimuTaskCreator creator )
{
    m_creators[ taskName ] = std::move( creator );
}

std::unique_ptr<ISimuTask> TaskRegistry::Create( const std::string& taskName ) const
{
    auto it = m_creators.find( taskName );
    if ( it == m_creators.end() )
    {
        return nullptr;
    }
    return it->second();
}

std::unique_ptr<ISimuTask> TaskRegistry::Create( TaskEnum task ) const
{
    return Create( TaskEnumToString( task ) );
}

bool TaskRegistry::Contains( const std::string& taskName ) const
{
    return m_creators.find( taskName ) != m_creators.end();
}

std::vector<std::string> TaskRegistry::GetAllRegisteredNames() const
{
    std::vector<std::string> names;
    names.reserve( m_creators.size() );
    for ( const auto& pair : m_creators )
    {
        names.push_back( pair.first );
    }
    return names;
}

const std::string& TaskEnumToString( TaskEnum task )
{
    static const std::string solve = "Solve";
    static const std::string grid = "Grid";
    static const std::string wallDist = "WallDist";
    static const std::string partition = "Partition";
    static const std::string theory = "Theory";
    static const std::string toyModel = "ToyModel";
    static const std::string postTask = "PostTask";
    static const std::string unknown = "Unknown";

    // Switch mapping tolerates retired numeric IDs without shifting the
    // remaining TaskEnum values.
    switch ( task )
    {
    case TaskEnum::SOLVE_FIELD:      return solve;
    case TaskEnum::CREATE_GRID:      return grid;
    case TaskEnum::CREATE_WALL_DIST: return wallDist;
    case TaskEnum::PARTITION_GRID:   return partition;
    case TaskEnum::SOLVE_THEORY:     return theory;
    case TaskEnum::TOY_MODEL:        return toyModel;
    case TaskEnum::POST_TASK:        return postTask;
    }
    return unknown;
}

EndNameSpace
