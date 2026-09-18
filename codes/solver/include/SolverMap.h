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
//#include "TextFileParser.h"
#include <map>
#include <memory>

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
    // Owning storage; GetSolver() returns non-owning Solver*.
    static HXVector< std::unique_ptr< Solver > > strSolver;
    static HXVector< std::unique_ptr< Solver > > unsSolver;
public:
    static void CreateSolvers();
    static void CreateSolvers( int gridType );
    // Injectable names (already U*/S* expanded). If null, uses SolverNameClass::GetSolverNames.
    // No numerical change: same BuildSolversInBucket path as the default overload.
    static void CreateSolvers( int gridType, const StringField * solverNameList );

    // S2 pure seam: injected non-null list wins; otherwise SolverNameClass for gridType.
    // No I/O, no SafeClone - unit-testable without full solver registry.
    static const StringField & SelectSolverNames(
        int gridType,
        const StringField * injected );

    static void FreeSolverMap();
    static void FreeSolverMap( int gridType );

    // Clear type↔index maps (and solverTypes). Called by FreeSolverMap;
    // also available for unit tests that exercise AddSolverInfo in isolation.
    static void ClearIndexMaps();

    static int GetSolverIndexBySolverType( int solverType );
    static int GetSolverTypeBySolverIndex( int solverIndex );
    static void AddSolverInfo( int solverType, int solverIndex );
    static Solver * GetSolver( int solverIndex, int gridType );
protected:
    static void AddSolverTypeToIndex( int solverType, int solverIndex );
    static void AddSolverIndexToType( int solverIndex, int solverType );
    static HXVector< std::unique_ptr< Solver > > * SolverBucket( int gridType );

    // S3: clone + StaticInit + index maps into the chosen bucket.
    // Does not touch SolverState / LusgsState (those stay in CreateSolvers).
    static void BuildSolversInBucket(
        int gridType,
        const StringField & solverNameList,
        HXVector< std::unique_ptr< Solver > > * solvers );
};

EndNameSpace
