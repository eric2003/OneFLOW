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

#include "Solver.h"
#include "SolverInfo.h"
#include "Fatal.h"
#include <map>
#include <memory>
#include <string>
#include <iostream>
#include <utility>


BeginNameSpace( ONEFLOW )

namespace
{
using SolverRegistry = std::map< std::string, std::unique_ptr< Solver > >;

SolverRegistry & GetSolverRegistry()
{
    static SolverRegistry registry;
    return registry;
}
}

Solver::Solver()
{
}

Solver::~Solver()
{
}

std::unique_ptr< Solver > Solver::SafeCloneUnique( const std::string & type )
{
    SolverRegistry & registry = GetSolverRegistry();
    SolverRegistry::iterator iter = registry.find( type );
    if ( iter == registry.end() )
    {
        Fatal( type + " class not found" );
        return nullptr;
    }

    return iter->second->Clone();
}


Solver * Solver::Register( const std::string & type, std::unique_ptr< Solver > clone )
{
    SolverRegistry & registry = GetSolverRegistry();
    SolverRegistry::iterator iter = registry.find( type );
    if ( iter != registry.end() ) return iter->second.get();

    Solver * registeredSolver = clone.get();
    registry.emplace( type, std::move( clone ) );
    return registeredSolver;
}

Solver * Solver::Register( const std::string & type, Solver * clone )
{
    return Solver::Register( type, std::unique_ptr< Solver >( clone ) );
}


EndNameSpace
