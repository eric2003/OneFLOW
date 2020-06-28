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

#include "AeroForceTaskReg.h"
#include "AeroForceTask.h"
#include "TaskState.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterForceTask )

void RegisterForceTask()
{
    REGISTER_DATA_CLASS( DumpAerodynamicForceTask );
    REGISTER_DATA_CLASS( DumpPressureCoeff );
}

void DumpAerodynamicForceTask( StringField & data )
{
    AerodynamicForceTask * task = new AerodynamicForceTask();
    TaskState::task = task;
}

void DumpPressureCoeff( StringField & data )
{
    int idump_pres = 1;
    CalcAeroForce( idump_pres );
}

EndNameSpace