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
#include "HXDefine.h"
#include "TextFileParser.h"
#include <map>

BeginNameSpace( ONEFLOW )

class Solver;
class SolverMap
{
public:
    SolverMap();
    ~SolverMap();
public:
    static IntField solverTypes;
    static std::map< int, int > solverTypeToIndex;
    static std::map< int, int > solverIndexToType;
    static HXVector< Solver * > strSolver;
    static HXVector< Solver * > unsSolver;
public:
    static void CreateSolvers();
    static void CreateSolvers( int gridType );
    static void FreeSolverMap();
    static void FreeSolverMap( int gridType );
    static int GetSolverIndexBySolverType( int solverType );
    static int GetSolverTypeBySolverIndex( int solverIndex );
    static void AddSolverInfo( int solverType, int solverIndex );
    static Solver * GetSolver( int solverIndex, int gridType );
protected:
    static void AddSolverTypeToIndex( int solverType, int solverIndex );
    static void AddSolverIndexToType( int solverIndex, int solverType );
};

class SolverNameClass
{
public:
    SolverNameClass();
    ~SolverNameClass();
public:
    static StringField unsSolverNameList;
    static StringField strSolverNameList;
    static bool flag;
public:
    static void Init();
    static void ReadSolverNames();
    static void ReadSolverNames( StringField & solverNameList );
    static StringField & GetSolverNames( int gridType );
};

EndNameSpace
