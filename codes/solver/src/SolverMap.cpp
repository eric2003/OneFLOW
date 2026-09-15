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

#include "SolverMap.h"
#include "Solver.h"
#include "TextFileParser.h"
#include "SolverNamePolicy.h"
#include "GridState.h"
#include "SolverState.h"
#include <map>
#include <iostream>


BeginNameSpace( ONEFLOW )

IntField SolverMap::solverTypes;
std::map< int, int > SolverMap::solverTypeToIndex;
std::map< int, int > SolverMap::solverIndexToType;
HXVector< Solver * > SolverMap::strSolver;
HXVector< Solver * > SolverMap::unsSolver;

SolverMap::SolverMap()
{
}

SolverMap::~SolverMap()
{
}

HXVector< Solver * > * SolverMap::SolverBucket( int gridType )
{
    if ( gridType == ONEFLOW::UMESH )
    {
        return & SolverMap::unsSolver;
    }
    return & SolverMap::strSolver;
}

void SolverMap::BuildSolversInBucket(
    int gridType,
    const StringField & solverNameList,
    HXVector< Solver * > * solvers )
{
    const int nSolver = static_cast< int >( solverNameList.size() );
    for ( int solverIndex = 0; solverIndex < nSolver; ++ solverIndex )
    {
        Solver * solver = Solver::SafeClone( solverNameList[ solverIndex ] );
        solver->solverIndex = solverIndex;
        solver->gridType = gridType;
        solver->StaticInit();

        SolverMap::AddSolverInfo( solver->solverType, solver->solverIndex );
        solvers->push_back( solver );
    }
}

void SolverMap::CreateSolvers( int gridType )
{
    // S1: select side (uns vs str)
    HXVector< Solver * > * solvers = SolverMap::SolverBucket( gridType );

    // S2: names from script/solver.txt + U/S prefix (SolverNamePolicy)
    StringField & solverNameList = SolverNameClass::GetSolverNames( gridType );
    const int nSolver = static_cast< int >( solverNameList.size() );

    LusgsState::Init( nSolver );

    // S3: SafeClone + StaticInit + index maps
    SolverMap::BuildSolversInBucket( gridType, solverNameList, solvers );

    // S4: solver-state side table
    SolverState::Init( nSolver );
}

void SolverMap::FreeSolverMap( int gridType )
{
    HXVector< Solver * > * solvers = SolverMap::SolverBucket( gridType );

    for ( int solverIndex = 0; solverIndex < static_cast< int >( solvers->size() ); ++ solverIndex )
    {
        delete ( * solvers )[ solverIndex ];
    }
    solvers->resize( 0 );
}

Solver * SolverMap::GetSolver( int solverIndex, int gridType )
{
    return ( * SolverMap::SolverBucket( gridType ) )[ solverIndex ];
}

void SolverMap::CreateSolvers()
{
    // S1: select side (uns vs str)
    // S2: names from script/solver.txt + U/S prefix
    // S3: SafeClone + StaticInit + index maps
    // S4: SolverState::Init
    // ... existing body unchanged ...
    SolverMap::CreateSolvers( ONEFLOW::UMESH );
    //SolverMap::CreateSolvers( ONEFLOW::SMESH );
}

void SolverMap::FreeSolverMap()
{
    SolverMap::FreeSolverMap( ONEFLOW::UMESH );
    SolverMap::FreeSolverMap( ONEFLOW::SMESH );
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

StringField SolverNameClass::unsSolverNameList;
StringField SolverNameClass::strSolverNameList;
bool SolverNameClass::flag = false;

SolverNameClass::SolverNameClass()
{
    ;
}

SolverNameClass::~SolverNameClass()
{
    ;
}

void SolverNameClass::Init()
{
    if ( flag ) return;
    flag = true;
    SolverNameClass::ReadSolverNames();
}

void SolverNameClass::ReadSolverNames()
{
    StringField solverNameList;
    SolverNameClass::ReadSolverNames( solverNameList );

    FillExpandedSolverNames(
        solverNameList,
        SolverNameClass::unsSolverNameList,
        SolverNameClass::strSolverNameList );
}

void SolverNameClass::ReadSolverNames( StringField & solverNameList )
{
    TextFileParser textFileParser;

    textFileParser.OpenPrjFile( "script/solver.txt", std::ios_base::in );

    // \t is the tab key
    const std::string keyWordSeparator = " ()\r\n\t#$,;\"";
    textFileParser.SetDefaultSeparator( keyWordSeparator );

    // Same pattern as MessageMapImp::ReadFile:
    // skip blank/comment lines; no spurious empty token at EOF.
    while ( textFileParser.ReadNextMeaningfulLine() )
    {
        std::string solverName = textFileParser.ReadNextWord();
        if ( solverName.empty() )
        {
            continue;
        }
        solverNameList.push_back( solverName );
    }

    textFileParser.CloseFile();
}

StringField & SolverNameClass::GetSolverNames( int gridType )
{
    SolverNameClass::Init();

    if ( gridType == ONEFLOW::UMESH )
    {
        return SolverNameClass::unsSolverNameList;
    }
    else
    {
        return SolverNameClass::strSolverNameList;
    }
}

void SolverNameClass::LoadFromBaseNames( const StringField & baseNames )
{
    SolverNameClass::unsSolverNameList.clear();
    SolverNameClass::strSolverNameList.clear();

    FillExpandedSolverNames(
        baseNames,
        SolverNameClass::unsSolverNameList,
        SolverNameClass::strSolverNameList );

    flag = true; // skip file path on later Init()
}

void SolverNameClass::Reset()
{
    unsSolverNameList.clear();
    strSolverNameList.clear();
    flag = false;
}


EndNameSpace
