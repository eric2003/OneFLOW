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
    // Index by underlying integer; keeps C++14/17 compilers happy without
    // a custom std::hash for enum class.
    static const std::string table[] = {
        "Solve",        // SOLVE_FIELD      = 0
        "Grid",         // CREATE_GRID      = 1
        "WallDist",     // CREATE_WALL_DIST = 2
        "Partition",    // PARTITION_GRID   = 3
        "FunctionTest", // FUNCTION_TEST    = 4
        "Theory",       // SOLVE_THEORY     = 5
        "ToyModel",     // TOY_MODEL        = 6
        "PostTask"      // POST_TASK        = 7
    };
    static const std::string unknown = "Unknown";

    const auto idx = static_cast<int>( task );
    if ( idx < 0 || idx >= static_cast<int>( sizeof( table ) / sizeof( table[ 0 ] ) ) )
    {
        return unknown;
    }
    return table[ idx ];
}

EndNameSpace
