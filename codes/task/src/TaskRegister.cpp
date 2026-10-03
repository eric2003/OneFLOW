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

#include "TaskRegister.h"
#include <memory>
#include <iostream>


BeginNameSpace( ONEFLOW )

std::unique_ptr< HXVector< VoidFunc > > TaskRegister::taskList;
std::unique_ptr< HXVector< std::string > > TaskRegister::taskNameList;

TaskRegister::TaskRegister()
{
}

TaskRegister::~TaskRegister()
{
}

void TaskRegister::Free()
{
    // reset is idempotent; subsequent Register() will re-allocate.
    TaskRegister::taskList.reset();
    TaskRegister::taskNameList.reset();
}

void TaskRegister::Register( VoidFunc taskfun, std::string const & taskname )
{
    if ( ! TaskRegister::taskList )
    {
        TaskRegister::taskList = std::make_unique< HXVector< VoidFunc > >();
        TaskRegister::taskNameList = std::make_unique< HXVector< std::string > >();
    }
    TaskRegister::taskList->push_back( taskfun );
    TaskRegister::taskNameList->push_back( taskname );
    //std::cout << "TaskRegister::Register " << taskname << "\n";
}

void TaskRegister::Run()
{
    // FIX: guard against Run() being called before any Register() call
    // has happened (taskList would still be null in that case).
    if ( ! TaskRegister::taskList )
    {
        return;
    }

    int n = TaskRegister::taskList->size();
    for ( int i = 0; i < n; ++ i )
    {
        VoidFunc & fun = ( * TaskRegister::taskList )[ i ];
        ( fun )( );
    }
}

class Tmp_Free_TaskRegister
{
public:
    Tmp_Free_TaskRegister() {}
    ~Tmp_Free_TaskRegister()
    {
        TaskRegister::Free();
    }
};

Tmp_Free_TaskRegister tmp_Free_TaskRegister;


EndNameSpace
