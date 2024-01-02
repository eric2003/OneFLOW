/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include <map>
#include <string>
#include <iostream>


BeginNameSpace( ONEFLOW )

std::map< std::string, Solver * > * Solver::classMap = 0;
Solver::Solver()
{
}

Solver::~Solver()
{
}

Solver * Solver::SafeClone( const std::string & type )
{
    std::map < std::string, Solver * >::iterator iter = Solver::classMap->find( type );
    if ( iter == Solver::classMap->end() )
    {
        std::cout << type << " class not found \n";
        exit( 0 );
    }

    return iter->second->Clone();
}

Solver * Solver::Register( const std::string & type, Solver * clone )
{
    if ( ! Solver::classMap )
    {
        Solver::classMap = new std::map < std::string, Solver * >();
    }

    std::map < std::string, Solver * >::iterator iter = Solver::classMap->find( type );
    if ( iter == Solver::classMap->end() )
    {
        ( * Solver::classMap )[ type ] = clone;
        return clone;
    }
    else
    {
        delete clone;
        return iter->second;
    }
}


EndNameSpace
