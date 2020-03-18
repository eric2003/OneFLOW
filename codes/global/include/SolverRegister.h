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
#include "SolverRegData.h"

BeginNameSpace( ONEFLOW )

class SolverRegData;
typedef SolverRegData * ( * SolverRegFun )( void );

#define REGISTER_REG_DATA( FUN ) \
class Init_SolverRegDataRegister##FUN \
{ \
public: \
    Init_SolverRegDataRegister##FUN() \
    {  \
        SolverRegister::Register( FUN ); \
    } \
};  \
Init_SolverRegDataRegister##FUN init_SolverRegDataRegister##FUN;

class SolverRegister
{
public:
    SolverRegister();
    ~SolverRegister();
public:
    static HXVector< SolverRegFun > * solverRegFunList;
public:
    static void Register( SolverRegFun solverRegFun );
    static void Run();
};

EndNameSpace