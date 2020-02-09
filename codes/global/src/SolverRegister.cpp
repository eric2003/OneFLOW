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
#include "SolverRegister.h"
#include "SolverRegData.h"
#include "SolverTaskReg.h"

BeginNameSpace( ONEFLOW )

HXVector< SolverRegFun > * SolverRegister::solverRegFunList = 0;

SolverRegister::SolverRegister()
{
}

SolverRegister::~SolverRegister()
{
}

void SolverRegister::Register( SolverRegFun solverRegFun )
{
    if ( ! SolverRegister::solverRegFunList )
    {
        SolverRegister::solverRegFunList = new HXVector< SolverRegFun >;
    }
    SolverRegister::solverRegFunList->push_back( solverRegFun );
}

void SolverRegister::Run()
{
    int n = SolverRegister::solverRegFunList->size();
    for ( int i = 0; i < n; ++ i )
    {
        SolverRegFun solverRegFun = ( * SolverRegister::solverRegFunList )[ i ];
        RegisterSolverTask( solverRegFun() );
    }
}

EndNameSpace