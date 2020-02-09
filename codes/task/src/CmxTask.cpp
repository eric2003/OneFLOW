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

#include "CmxTask.h"
#include "SimpleTask.h"
#include "Command.h"
#include "Register.h"
#include "Message.h"
#include "SolverState.h"
#include "SolverDef.h"
#include "GridState.h"
#include "HXClone.h"
#include "Task.h"
#include "TaskState.h"
#include "Solver.h"
#include "Category.h"
#include "SolverMap.h"
#include "BgGrid.h"
#include "Zone.h"
#include "Grid.h"
#include "LogFile.h"

BeginNameSpace( ONEFLOW )

void CmdBasicAction( int funcType )
{
    SolverState::msgId = TaskState::task->taskId;
    HXClone * cloneClass = ONEFLOW::GetClass( SolverState::msgId, SolverState::tid, funcType );
    if ( cloneClass )
    {
        cloneClass->Solve();
    }
}

void CmdAction()
{
    CmdBasicAction( COMM_FUNC );
}

void CmdActionNext()
{
    CmdBasicAction( RECV_FUNC );
}

void GenerateCmdList( int msgId )
{
    int sTid = SolverState::tid;

    HXRegister * hxRegister = RegisterFactory::GetRegister( sTid, MESG_FUNC );

    string msgName = MessageMap::GetMsgName( msgId );

    HXClone * cloneClass = hxRegister->GetClass( msgName );

    if ( cloneClass )
    {
        cloneClass->Solve();
    }
    else
    {
        ONEFLOW::AddCmdToList( msgName );
    }
}

void AddCmdToList( const string & msgName )
{
    int msgId = MessageMap::GetMsgId( msgName );

    ONEFLOW::AddCmdToList( msgId, SolverState::tid );
}

void AddCmdToList( int msgId, int sTid )
{
    ONEFLOW::CreateTask( msgId, sTid );

    ONEFLOW::SetFile( msgId, sTid );

    SimpleCmd * cmd = new SimpleCmd();

    cmd->AddTask( TaskState::task );

    CMD::AddCmd( cmd );
}

void SetTaskAction()
{
    TaskState::task->action     = & ONEFLOW::CmdAction;
    TaskState::task->sendAction = & ONEFLOW::CmdAction;
    TaskState::task->recvAction = & ONEFLOW::CmdActionNext;
}

void CreateTask( int msgId, int sTid )
{
    HXClone * cloneClass = ONEFLOW::GetClass( msgId, sTid, TASK_FUNC );

    if ( cloneClass )
    {
        cloneClass->Solve();
    }
    else
    {
        TaskState::task = new SimpleTask();
    }

    TaskState::task->taskId = msgId;
    TaskState::task->taskName = MessageMap::GetMsgName( msgId );

    SolverState::tid = sTid;
    SetTaskAction();
}

void SetFile( int msgId, int sTid )
{
    HXClone * cloneClass = ONEFLOW::GetClass( msgId, sTid, FILE_FUNC );

    if ( cloneClass )
    {
        cloneClass->Solve();
    }
}

HXClone * GetClass( int msgId, int sTid, int msgType )
{
    HXRegister * hxRegister = RegisterFactory::GetRegister( sTid, msgType );

    string msgName = MessageMap::GetMsgName( msgId );

    HXClone * cloneClass = hxRegister->GetClass( msgName );

    return cloneClass;
}

void SsSgTask( const string & taskName )
{
    int taskCode = MessageMap::GetMsgId( taskName );

    ONEFLOW::GenerateCmdList( taskCode );

    CMD::ExecuteCmd();
}

void MsMgTask( const string & taskname )
{
    for ( int sid = 0; sid < SolverState::nSolver; ++ sid )
    {
        SolverState::SetTidById( sid );

        for ( int gl = 0; gl < GridState::nGrids; ++ gl )
        {
            GridState::SetGridLevel( gl );
            ONEFLOW::SsSgTask( taskname );
        }
    }
}


EndNameSpace