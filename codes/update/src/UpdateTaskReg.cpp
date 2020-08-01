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

#include "UpdateTaskReg.h"
#include "Update.h"
#include "SolverState.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterUpdateTask )

REGISTER_TASK(RegisterINsUpdateTask)


void UpdateFlowField( StringField & data )
{
    int sTid = SolverState::tid;
    Update * update = CreateUpdate( sTid );
    update->UpdateFlowField( sTid );
    delete update;
}

void UpdateINsFlowField(StringField & data)
{
	int sTid = SolverState::tid;
	Update * update = CreateUpdate(sTid);
	update->UpdateINsFlowField(sTid);
	delete update;
}

void RegisterUpdateTask()
{
	REGISTER_DATA_CLASS(UpdateFlowField);
}

void RegisterINsUpdateTask()
{
	REGISTER_DATA_CLASS(UpdateINsFlowField);
}



EndNameSpace


