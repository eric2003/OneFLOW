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

#include "TaskImp.h"
#include "TaskCom.h"
#include "TaskState.h"
#include "ReadTask.h"
#include "WriteTask.h"
#include "InterfaceTask.h"
#include "OversetTask.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterComTask )

void RegisterComTask()
{
    REGISTER_DATA_CLASS( ReadBinaryFileTask  );
    REGISTER_DATA_CLASS( WriteBinaryFileTask );
    REGISTER_DATA_CLASS( ServerUpdateInterfaceTask  );
    REGISTER_DATA_CLASS( WriteAsciiFileTask );
    REGISTER_DATA_CLASS( ServerUpdateOversetInterfaceTask );
}

void ReadBinaryFileTask( StringField & data )
{
    CReadFile * task = new CReadFile();
    task->mainAction = & ReadBinaryFile;
    TaskState::task = task;
}

void WriteBinaryFileTask( StringField & data )
{
    CWriteFile * task = new CWriteFile();
    task->mainAction = & WriteBinaryFile;
    TaskState::task = task;
}

void WriteAsciiFileTask( StringField & data )
{
    CWriteFile * task = new CWriteFile();
    task->mainAction = & WriteAsciiFile;
    TaskState::task = task;
}

void ServerUpdateInterfaceTask( StringField & data )
{
    CUpdateInterface * task = new CUpdateInterface();
    TaskState::task = task;
}

void ServerUpdateOversetInterfaceTask( StringField & data )
{
    OversetTask * task = new OversetTask();
    TaskState::task = task;
}

EndNameSpace