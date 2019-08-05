/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "RegData.h"
#include "SolverTaskReg.h"

BeginNameSpace( ONEFLOW )

HXVector< RegDataFun > * RegDataRegister::regDataFunList = 0;

RegDataRegister::RegDataRegister()
{
}

RegDataRegister::~RegDataRegister()
{
}

void RegDataRegister::Register( RegDataFun regDataFun )
{
    if ( ! RegDataRegister::regDataFunList )
    {
        RegDataRegister::regDataFunList = new HXVector< RegDataFun >;
    }
    RegDataRegister::regDataFunList->push_back( regDataFun );
}

void RegDataRegister::Run()
{
    int n = RegDataRegister::regDataFunList->size();
    for ( int i = 0; i < n; ++ i )
    {
        RegDataFun regDataFun = ( * RegDataRegister::regDataFunList )[ i ];
        RegisterSolverTask( regDataFun() );
    }
}

EndNameSpace