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
// Index maps + SelectSolverNames only (no SafeClone, no SolverState, no Lusgs).
// Linked by production SolverMap.cpp consumers and by solver_map_index_test.

#include "SolverMap.h"
#include "SolverNameList.h"

BeginNameSpace( ONEFLOW )

IntField SolverMap::solverTypes;
std::map< int, int > SolverMap::solverTypeToIndex;
std::map< int, int > SolverMap::solverIndexToType;

void SolverMap::ClearIndexMaps()
{
    SolverMap::solverTypes.clear();
    SolverMap::solverTypeToIndex.clear();
    SolverMap::solverIndexToType.clear();
}

const StringField & SolverMap::SelectSolverNames(
    int gridType,
    const StringField * injected )
{
    if ( injected != nullptr )
    {
        return *injected;
    }
    return SolverNameClass::GetSolverNames( gridType );
}

int SolverMap::GetSolverIndexBySolverType( int solverType )
{
    std::map< int, int >::iterator iter;
    iter = SolverMap::solverTypeToIndex.find( solverType );
    return iter->second;
}

int SolverMap::GetSolverTypeBySolverIndex( int solverIndex )
{
    std::map< int, int >::iterator iter;
    iter = SolverMap::solverIndexToType.find( solverIndex );
    return iter->second;
}

void SolverMap::AddSolverInfo( int solverType, int solverIndex )
{
    SolverMap::AddSolverTypeToIndex( solverType, solverIndex );
    SolverMap::AddSolverIndexToType( solverIndex, solverType );
}

void SolverMap::AddSolverTypeToIndex( int solverType, int solverIndex )
{
    std::map< int, int >::iterator iter;
    iter = SolverMap::solverTypeToIndex.find( solverType );
    if ( iter == SolverMap::solverTypeToIndex.end() )
    {
        SolverMap::solverTypeToIndex[ solverType ] = solverIndex;
        SolverMap::solverTypes.push_back( solverType );
    }
}

void SolverMap::AddSolverIndexToType( int solverIndex, int solverType )
{
    std::map< int, int >::iterator iter = SolverMap::solverIndexToType.find( solverIndex );
    if ( iter == SolverMap::solverIndexToType.end() )
    {
        SolverMap::solverIndexToType[ solverIndex ] = solverType;
    }
}

EndNameSpace
