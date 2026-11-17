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

#include "TaskImp.h"
#include <memory>
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
    auto task = std::make_unique<CReadFile>();
    task->mainAction = & ReadBinaryFile;
    TaskState::createdTask = std::move( task );
}

void WriteBinaryFileTask( StringField & data )
{
    auto task = std::make_unique<CWriteFile>();
    task->mainAction = & WriteBinaryFile;
    TaskState::createdTask = std::move( task );
}

void WriteAsciiFileTask( StringField & data )
{
    auto task = std::make_unique<CWriteFile>();
    task->mainAction = & WriteAsciiFile;
    TaskState::createdTask = std::move( task );
}

void ServerUpdateInterfaceTask( StringField & data )
{
    TaskState::createdTask = std::make_unique<CUpdateInterface>();
}

void ServerUpdateOversetInterfaceTask( StringField & data )
{
    TaskState::createdTask = std::make_unique<OversetTask>();
}

EndNameSpace
