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


#pragma once
#include "HXDefine.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

class SolverP;
class Solver;
class LusgsSolver;
class LusgsPair;

class LusgsState
{
public:
    LusgsState();
    ~LusgsState();
public:
    static HXVector< LusgsSolver * > str;
    static HXVector< LusgsSolver * > uns;
public:
    static void Init( int nSolver );
    static void AddSolver( int sid, int gridType, LusgsSolver * solver );
    static LusgsSolver * GetLusgsSolver();
};

class SolverState
{
public:
    SolverState();
    ~SolverState();
public:
    static int id;
    static int tid;
    static int nSolver;
    static int msgId;
    static IntField convergeFlag;
public:
    static void Init( int nSolver );
    static void SetTid( int tid );
    static void SetTidById( int id );
    static Solver * GetSolver();
public:
    static bool Converge();
};


EndNameSpace