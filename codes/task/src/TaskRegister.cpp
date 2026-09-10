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
#include <iostream>


BeginNameSpace( ONEFLOW )

HXVector< VoidFunc > * TaskRegister::taskList = 0;
HXVector< std::string > * TaskRegister::taskNameList = 0;

TaskRegister::TaskRegister()
{
}

TaskRegister::~TaskRegister()
{
}

//void TaskRegister::Free()
//{
//    delete TaskRegister::taskList;
//    delete TaskRegister::taskNameList;
//}

void TaskRegister::Free()
{
    delete TaskRegister::taskList;
    delete TaskRegister::taskNameList;

    // FIX: null out pointers after delete. Without this, Free() leaves
    // dangling pointers, which causes two separate hazards:
    //   1. A subsequent Register() call would push_back into freed
    //      memory (use-after-free), since the null-check in Register()
    //      only guards against a NULL pointer, not a dangling one.
    //   2. Free() itself was not idempotent: calling it twice (e.g. once
    //      manually, then again via Tmp_Free_TaskRegister's destructor
    //      at program exit) would double-free the same memory.
    // Nulling the pointers makes Free() safe to call multiple times and
    // makes the class safely re-usable after being freed.
    TaskRegister::taskList = 0;
    TaskRegister::taskNameList = 0;
}

void TaskRegister::Register( VoidFunc taskfun, std::string const & taskname )
{
    if ( ! TaskRegister::taskList )
    {
        TaskRegister::taskList = new HXVector< VoidFunc >;
        TaskRegister::taskNameList = new HXVector< std::string >;
    }
    TaskRegister::taskList->push_back( taskfun );
    TaskRegister::taskNameList->push_back( taskname );
    //std::cout << "TaskRegister::Register " << taskname << "\n";
}

//void TaskRegister::Run()
//{
//    int n = TaskRegister::taskList->size();
//    for ( int i = 0; i < n; ++ i )
//    {
//        VoidFunc & fun = ( * TaskRegister::taskList )[ i ];
//        ( fun )( );
//    }
//}

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
