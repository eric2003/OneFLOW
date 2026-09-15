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
#include "HXDefine.h"   // if StringField is used
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

// script/solver.txt base name (e.g. "NsSolver");
// prefix by grid type, then Solver::SafeClone

inline std::string MakeUnstructuredSolverName( const std::string& baseName )
{
    return "U" + baseName;
}

inline std::string MakeStructuredSolverName( const std::string& baseName )
{
    return "S" + baseName;
}

inline std::vector<std::string> ExpandSolverNamesForUnstructured(
    const std::vector<std::string>& baseNames )
{
    std::vector<std::string> out;
    out.reserve( baseNames.size() );
    for ( const auto& base : baseNames )
    {
        out.push_back( MakeUnstructuredSolverName( base ) );
    }
    return out;
}

inline std::vector<std::string> ExpandSolverNamesForStructured(
    const std::vector<std::string>& baseNames )
{
    std::vector<std::string> out;
    out.reserve( baseNames.size() );
    for ( const auto& base : baseNames )
    {
        out.push_back( MakeStructuredSolverName( base ) );
    }
    return out;
}

// Must come AFTER Make* ¡ª GCC two-phase lookup requires visible declarations
// at the point of the template definition (MSVC is more permissive).
template< typename StringRange >
inline void FillExpandedSolverNames(
    const StringRange& baseNames,
    StringField& outUns,
    StringField& outStr )
{
    outUns.clear();
    outStr.clear();
    outUns.reserve( baseNames.size() );
    outStr.reserve( baseNames.size() );
    for ( const auto& base : baseNames )
    {
        outUns.push_back( MakeUnstructuredSolverName( base ) );
        outStr.push_back( MakeStructuredSolverName( base ) );
    }
}

EndNameSpace