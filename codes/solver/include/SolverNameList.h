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


BeginNameSpace( ONEFLOW )

// Holds expanded solver registration names (U*/S*) for UMESH / SMESH.
// File I/O + SolverNamePolicy expansion; no Solver* / SafeClone.
class SolverNameClass
{
public:
    SolverNameClass();
    ~SolverNameClass();

    static StringField unsSolverNameList;
    static StringField strSolverNameList;
    static bool flag;

    static void Init();
    static void ReadSolverNames();
    static void ReadSolverNames( StringField & solverNameList );
    static StringField & GetSolverNames( int gridType );

    // Test / injection: expand base names, no script/solver.txt
    static void LoadFromBaseNames( const StringField & baseNames );
    static void Reset();
};

EndNameSpace
