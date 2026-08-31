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

IntField SolverMap::tid;
std::map< int, int > SolverMap::tid2Id;
std::map< int, int > SolverMap::id2Tid;
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

int SolverMap::GetId( int solverType )
{
    std::map< int, int >::iterator iter;
    iter = SolverMap::tid2Id.find( solverType );
    return iter->second;
}

int SolverMap::GetTid( int solverIndex )
{
    std::map< int, int >::iterator iter;
    iter = SolverMap::id2Tid.find( solverIndex );
    return iter->second;
}

void SolverMap::AddSolverInfo( int solverType, int solverIndex )
{
    SolverMap::AddTid2Id( solverType, solverIndex );
    SolverMap::AddId2Tid( solverIndex, solverType );
}

void SolverMap::AddTid2Id( int solverType, int solverIndex )
{
    std::map< int, int >::iterator iter;
    iter = SolverMap::tid2Id.find( solverType );
    if ( iter == SolverMap::tid2Id.end() )
    {
        SolverMap::tid2Id[ solverType ] = solverIndex;
        SolverMap::tid.push_back( solverType );
    }
}

void SolverMap::AddId2Tid( int solverIndex, int solverType )
{
    std::map< int, int >::iterator iter = SolverMap::id2Tid.find( solverIndex );
    if ( iter == SolverMap::id2Tid.end() )
    {
        SolverMap::id2Tid[ solverIndex ] = solverType;
    }
}

Solver * SolverMap::GetSolver( int id, int gridType )
{
    if ( gridType == ONEFLOW::UMESH )
    {
        return unsSolver[ id ];
    }
    else
    {
        return strSolver[ id ];
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
