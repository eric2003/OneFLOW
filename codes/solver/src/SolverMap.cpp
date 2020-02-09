/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
using namespace std;

BeginNameSpace( ONEFLOW )

IntField SolverMap::tid;
map< int, int > SolverMap::tid2Id;
map< int, int > SolverMap::id2Tid;
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
    for ( int sid = 0; sid < nSolver; ++ sid )
    {
        Solver * solver = Solver::SafeClone( solverNameList[ sid ] );
        solver->sid = sid;
        solver->gridType = gridType;
        solver->StaticInit();
        
        SolverMap::AddSolverInfo( solver->sTid, solver->sid );
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

    for ( int sid = 0; sid < solvers->size(); ++ sid )
    {
        Solver * solver = ( * solvers )[ sid ];
        delete solver;
    }
    solvers->resize( 0 );
}

void SolverMap::FreeSolverMap()
{
    SolverMap::FreeSolverMap( ONEFLOW::UMESH );
    SolverMap::FreeSolverMap( ONEFLOW::SMESH );
}

int SolverMap::GetId( int sTid )
{
    map< int, int >::iterator iter;
    iter = SolverMap::tid2Id.find( sTid );
    return iter->second;
}

int SolverMap::GetTid( int sid )
{
    map< int, int >::iterator iter;
    iter = SolverMap::id2Tid.find( sid );
    return iter->second;
}

void SolverMap::AddSolverInfo( int sTid, int sid )
{
    SolverMap::AddTid2Id( sTid, sid );
    SolverMap::AddId2Tid( sid, sTid );
}

void SolverMap::AddTid2Id( int sTid, int sid )
{
    map< int, int >::iterator iter;
    iter = SolverMap::tid2Id.find( sTid );
    if ( iter == SolverMap::tid2Id.end() )
    {
        SolverMap::tid2Id[ sTid ] = sid;
        SolverMap::tid.push_back( sTid );
    }
}

void SolverMap::AddId2Tid( int sid, int sTid )
{
    map< int, int >::iterator iter = SolverMap::id2Tid.find( sid );
    if ( iter == SolverMap::id2Tid.end() )
    {
        SolverMap::id2Tid[ sid ] = sTid;
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
        string solverName = solverNameList[ isol ];

        ONEFLOW::StrIO.ClearAll();
        ONEFLOW::StrIO << "U" << solverName;
        string uSolverName = ONEFLOW::StrIO.str();

        ONEFLOW::StrIO.ClearAll();
        ONEFLOW::StrIO << "S" << solverName;
        string sSolverName = ONEFLOW::StrIO.str();

        SolverNameClass::unsSolverNameList.push_back( uSolverName );
        SolverNameClass::strSolverNameList.push_back( sSolverName );
    }
}

void SolverNameClass::ReadSolverNames( StringField & solverNameList )
{
    FileIO ioFile;

    ioFile.OpenPrjFile( "script/solver.txt", ios_base::in );

    //\tÎªtab¼ü
    string keyWordSeparator = " ()\r\n\t#$,;\"";
    ioFile.SetDefaultSeparator( keyWordSeparator );

    while ( ! ioFile.ReachTheEndOfFile()  )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        string solverName = ioFile.ReadNextWord();
        solverNameList.push_back( solverName );
    }

    ioFile.CloseFile();
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