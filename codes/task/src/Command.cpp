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

#include "Command.h"
#include "Task.h"
#include "TaskState.h"
#include "TimeSpan.h"
#include <iostream>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

Command::Command()
{
    tasks = new TList;
}

Command::~Command()
{
    for ( UInt iTask = 0; iTask < tasks->size(); ++ iTask )
    {
        delete ( * tasks )[ iTask ];
    }
    delete tasks;
}

void Command::AddTask( Task * task )
{
    tasks->push_back( task );
}

SimpleCmd::SimpleCmd()
{
}

SimpleCmd::~SimpleCmd()
{
}

void SimpleCmd::Execute()
{
    for ( UInt iTask = 0; iTask < tasks->size(); ++ iTask )
    {
        Task * task = ( * tasks )[ iTask ];
        TaskState::task = task;
        task->Run();
    }
}

HXVector< Command * > * CMD::cmdList = 0;

CMD::CMD()
{
    ;
}

CMD::~CMD()
{
    ;
}

void CMD::Init()
{
    if ( CMD::cmdList ) return;
    cmdList = new HXVector< Command * >;
}

void CMD::Free()
{
    delete cmdList;
    cmdList = 0;
}

void CMD::AddCmd( Command * cmd )
{
    CMD::Init();
    CMD::cmdList->push_back( cmd );
}

void CMD::RunCmd( Command * cmd )
{
    cmd->Execute();
}

void CMD::Clear()
{
    CMD::cmdList->resize( 0 );
}

void CMD::ExecuteCmd()
{
    UInt nCmd = CMD::cmdList->size();
    for ( UInt iCmd = 0; iCmd < nCmd; ++ iCmd )
    {
        Command * cmd = ( * CMD::cmdList )[ iCmd ];
        //CMD::ShowCmdInfo( cmd, iCmd );
        cmd->Execute();
        delete cmd;
    }
    CMD::Clear();
}

void CMD::ShowCmdInfo( Command * cmd, int iCmd )
{
    HXVector< Task * > * tasks = cmd->tasks;
    for ( UInt i = 0; i < tasks->size(); ++ i )
    {
        Task * task = ( * tasks ) [ i ];
        int iTaskGlobal = iCmd + i;

        cout << " iTaskGlobal = " << iTaskGlobal << " iTaskLocal = " << i << " ";
        cout << " TaskCode = " << task->taskId << " Task Name = " << task->taskName << endl;
    }
}


EndNameSpace