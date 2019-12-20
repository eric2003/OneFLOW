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


#pragma once
#include "HXDefine.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#define REGISTER_TASK( TASK ) \
class Init_Register##TASK \
{ \
public: \
    Init_Register##TASK() \
    {  \
        TaskRegister::Register( TASK ); \
    } \
};  \
Init_Register##TASK init_Register##TASK;

class TaskRegister
{
public:
    TaskRegister();
    ~TaskRegister();
public:
    static HXVector< VoidFunc > * taskList;
public:
    static void Register( VoidFunc taskfun );
    static void Run();
};

EndNameSpace