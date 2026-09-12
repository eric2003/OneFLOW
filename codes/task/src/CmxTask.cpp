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
#include "Zone.h"
#include "Grid.h"
#include "LogFile.h"
#include <memory>
#include <utility>

BeginNameSpace( ONEFLOW )

HXClone * GetClass(
    int operationId,
    int solverType,
    int funcType )
{
    HXRegister * hxRegister =
        RegisterFactory::GetRegister(
            solverType,
            funcType );

    const std::string operationName =
        MessageMap::GetMsgName( operationId );

    HXClone * cloneClass =
        hxRegister->GetClass( operationName );

    return cloneClass;
}

// ============================================================
// Operation planning
// ============================================================

void GenerateCmdList( int operationId )
{
    const int solverType = SolverState::solverType;

    HXRegister * hxRegister =
        RegisterFactory::GetRegister( solverType, MESG_FUNC );

    const std::string operationName =
        MessageMap::GetMsgName( operationId );

    HXClone * cloneClass =
        hxRegister->GetClass( operationName );

    if ( cloneClass )
    {
        // Expand the operation into runtime commands.
        cloneClass->Solve();
    }
    else
    {
        // Treat the operation as a single runtime command.
        ONEFLOW::AddCmdToList( operationName );
    }
}

// ============================================================
// Command construction
// ============================================================

void AddCmdToList( const std::string & operationName )
{
    const int operationId =
        MessageMap::GetMsgId( operationName );

    ONEFLOW::AddCmdToList(
        operationId,
        SolverState::solverType );
}

void AddCmdToList(
    int operationId,
    int solverType )
{
    // Build the task associated with the operation.
    Task * task =
        ONEFLOW::CreateTask(
            operationId,
            solverType );

    if ( task == nullptr )
    {
        return;
    }

    // Prepare files/resources required by the operation.
    ONEFLOW::SetFile(
        task,
        operationId,
        solverType );

    // Build the command with RAII ownership.
    std::unique_ptr< SimpleCmd > cmd(
        new SimpleCmd() );

    cmd->AddTask( task );

    // Transfer command ownership to CMD.
    CMD::AddCmd( std::move( cmd ) );
}

// ============================================================
// Task construction
// ============================================================

Task * CreateTask( int operationId, int solverType )
{
    HXClone * cloneClass =
        ONEFLOW::GetClass(
            operationId,
            solverType,
            TASK_FUNC );

    if ( cloneClass )
    {
        cloneClass->Solve();
    }
    else
    {
        // Use the default task implementation.
        TaskState::task = new SimpleTask();
    }

    Task * task = TaskState::task;

    task->taskId = operationId;
    TaskState::task->taskName =
        MessageMap::GetMsgName( operationId );

    SolverState::solverType = solverType;

    SetTaskAction( task );

    return task;
}


// ============================================================
// Resource preparation
// ============================================================

void SetFile( Task * task, int operationId, int solverType )
{
    HXClone * cloneClass =
        ONEFLOW::GetClass(
            operationId,
            solverType,
            FILE_FUNC );

    if ( cloneClass )
    {
        // Keep the current TaskState mechanism temporarily.
        cloneClass->Solve();
    }
}

// ============================================================
// Action dispatch
// ============================================================

void SetTaskAction(Task * task)
{
    if (task == nullptr)
    {
        return;
    }

    task->action = CmdAction;
    task->sendAction = CmdAction;
    task->recvAction = CmdActionNext;
}

void CmdBasicAction( int funcType )
{
    SolverState::msgId =
        TaskState::task->taskId;

    HXClone * cloneClass =
        ONEFLOW::GetClass(
            SolverState::msgId,
            SolverState::solverType,
            funcType );

    if ( cloneClass )
    {
        // Execute the concrete registered implementation.
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



// ============================================================
// Operation execution entry
// ============================================================

void SingleSolverSingleGridTask( const std::string & taskName )
{
    // Resolve the operation name.
    const int operationId =
        MessageMap::GetMsgId( taskName );

    // Build the execution plan for the operation.
    GenerateCmdList( operationId );

    // Execute the generated plan.
    CMD::ExecuteCmd();
}


// ============================================================
// Multi-solver / multi-grid execution
// ============================================================

void MultiSolverMultiGridTask( const std::string & taskName )
{
    for ( int solverIndex = 0;
        solverIndex < SolverState::nSolver;
        ++ solverIndex )
    {
        SolverState::SetSolverTypeBySolverIndex( solverIndex );

        for ( int gl = 0;
            gl < GridState::nGrids;
            ++ gl )
        {
            GridState::SetGridLevel( gl );

            ONEFLOW::SingleSolverSingleGridTask(
                taskName );
        }
    }
}


EndNameSpace