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
#include "GridState.h"
#include "SolverState.h"
#include "OStream.h"
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

void SolverMap::CreateSolvers()
{
    SolverMap::CreateSolvers( ONEFLOW::UMESH );
    //SolverMap::CreateSolvers( ONEFLOW::SMESH );
}

void SolverMap::CreateSolvers( int gridType )
{
    HXVector< Solver * > * solvers = 0;
    if ( gridType == ONEFLOW::UMESH )
    {
        solvers = & SolverMap::unsSolver;
    }
    else
    {
        solvers = & SolverMap::strSolver;
    }

    StringField & solverNameList = SolverNameClass::GetSolverNames( gridType );
    int nSolver = solverNameList.size();

    LusgsState::Init( nSolver );
    for ( int solverIndex = 0; solverIndex < nSolver; ++ solverIndex )
    {
        Solver * solver = Solver::SafeClone( solverNameList[ solverIndex ] );
        solver->solverIndex = solverIndex;
        solver->gridType = gridType;
        solver->StaticInit();
        
        SolverMap::AddSolverInfo( solver->solverType, solver->solverIndex );
        solvers->push_back( solver );
    }

    SolverState::Init( nSolver );
}

void SolverMap::FreeSolverMap( int gridType )
{
    HXVector< Solver * > * solvers = 0;
    if ( gridType == ONEFLOW::UMESH )
    {
        solvers = & SolverMap::unsSolver;
    }
    else
    {
        solvers = & SolverMap::strSolver;
    }

    for ( int solverIndex = 0; solverIndex < solvers->size(); ++ solverIndex )
    {
        Solver * solver = ( * solvers )[ solverIndex ];
        delete solver;
    }
    solvers->resize( 0 );
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

Solver * SolverMap::GetSolver( int solverIndex, int gridType )
{
    if ( gridType == ONEFLOW::UMESH )
    {
        return unsSolver[ solverIndex ];
    }
    else
    {
        return strSolver[ solverIndex ];
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
    for ( int isol = 0; isol < solverNameList.size(); ++ isol )
    {
        std::string solverName = solverNameList[ isol ];

        OStream &logger = OStream::Instance();
        logger.ClearAll();
        logger << "U" << solverName;
        std::string uSolverName = logger.str();

        logger.ClearAll();
        logger << "S" << solverName;
        std::string sSolverName = logger.str();

        SolverNameClass::unsSolverNameList.push_back( uSolverName );
        SolverNameClass::strSolverNameList.push_back( sSolverName );
    }
}

void SolverNameClass::ReadSolverNames( StringField & solverNameList )
{
    TextFileParser textFileParser;

    textFileParser.OpenPrjFile( "script/solver.txt", std::ios_base::in );

    //\t is the tab key
    std::string keyWordSeparator = " ()\r\n\t#$,;\"";
    textFileParser.SetDefaultSeparator( keyWordSeparator );

    while ( ! textFileParser.ReachTheEndOfFile()  )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        std::string solverName = textFileParser.ReadNextWord();
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


EndNameSpace
