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
#include <map>
using namespace std;
BeginNameSpace( ONEFLOW )

class Solver;
#define IMPLEMENT_SOLVER_CLONE( TYPE ) \
Solver * Clone() const { return new TYPE( * this ); }

#define REGISTER_SOLVER( TYPE ) \
    Solver * TYPE ## _myClass = \
        Solver::Register( #TYPE, new TYPE() );

class SolverInfo;

class Solver
{
public:
    Solver();
    virtual ~Solver();
public:
    virtual Solver * Clone() const = 0;
public:
    static Solver * SafeClone( const string & type );
    static Solver * Register( const string & type, Solver * clone );
    static map < string, Solver * > * classMap;
public:
    int sid, sTid;
    int gridType;
public:
    virtual void StaticInit(){};
};

EndNameSpace