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

#include "TaskRegister.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

HXVector< VoidFunc > * TaskRegister::taskList = 0;

TaskRegister::TaskRegister()
{
}

TaskRegister::~TaskRegister()
{
}

void TaskRegister::Register( VoidFunc taskfun )
{
    if ( ! TaskRegister::taskList )
    {
        TaskRegister::taskList = new HXVector< VoidFunc >;
    }
    TaskRegister::taskList->push_back( taskfun );
}

void TaskRegister::Run()
{
    int n = TaskRegister::taskList->size();
    for ( int i = 0; i < n; ++ i )
    {
        VoidFunc & fun = ( * TaskRegister::taskList )[ i ];
        ( fun )( );
    }
}


EndNameSpace