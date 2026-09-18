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
#include "FileMap.h"
#include <memory>
#include <utility>
#include <map>
#include <cstdint>

BeginNameSpace( ONEFLOW )

namespace {

// Cache HXClone* by (operationId, solverType, funcType). Invalidated when
// MessageMap::Epoch() changes (Init/Free). Hot path: CmdAction / GenerateCmdList.
struct GetClassCache
{
    int epoch = -1;
    std::map< std::uint64_t, HXClone * > table;

    static std::uint64_t MakeKey( int operationId, int solverType, int funcType )
    {
        return ( static_cast< std::uint64_t >( static_cast< std::uint32_t >( operationId ) ) )
             | ( static_cast< std::uint64_t >( static_cast< std::uint32_t >( solverType ) ) << 20 )
             | ( static_cast< std::uint64_t >( static_cast< std::uint32_t >( funcType ) ) << 40 );
    }

    void SyncEpoch()
    {
        const int ep = MessageMap::Epoch();
        if ( ep != epoch )
        {
            table.clear();
            epoch = ep;
        }
    }
};

GetClassCache & ClassCache()
{
    static GetClassCache cache;
    return cache;
}

} // namespace

HXClone * GetClass(
    int operationId,
    int solverType,
    int funcType )
{
    GetClassCache & cache = ClassCache();
    cache.SyncEpoch();

    const std::uint64_t key =
        GetClassCache::MakeKey( operationId, solverType, funcType );

    auto it = cache.table.find( key );
    if ( it != cache.table.end() )
    {
        return it->second; // may be nullptr (negative cache)
    }

    HXRegister * hxRegister =
        RegisterFactory::GetRegister(
            solverType,
            funcType );

    const std::string & operationName =
        MessageMap::GetMsgName( operationId );

    HXClone * cloneClass =
        hxRegister->GetClass( operationName );

    cache.table[ key ] = cloneClass;
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

    const std::string & operationName =
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
    ONEFLOW::ConfigureTaskFile(
        task,
        operationId,
        solverType );

    // Take temporary ownership of the newly created Task.
    std::unique_ptr< Task > ownedTask( task );

    // Build the command with RAII ownership.
    std::unique_ptr< SimpleCmd > cmd(
        new SimpleCmd() );

    // Transfer Task ownership to the Command.
    cmd->AddTask( std::move( ownedTask ) );

    // Transfer Command ownership to CMD.
    CMD::AddCmd( std::move( cmd ) );
}


// ============================================================
// Task construction
// ============================================================

namespace
{

    Task * CreateTaskByRegisteredFunction(
        HXClone * cloneClass )
    {
        if ( cloneClass == nullptr )
        {
            return nullptr;
        }

        // TASK_FUNC callbacks return their construction result here.
        TaskState::createdTask = nullptr;

        cloneClass->Solve();

        Task * task = TaskState::createdTask;

        // Do not keep a stale construction result.
        TaskState::createdTask = nullptr;

        return task;
    }

}

Task * CreateTask( int operationId, int solverType )
{
    Task * task = nullptr;

    HXClone * cloneClass =
        ONEFLOW::GetClass(
            operationId,
            solverType,
            TASK_FUNC );

    if ( cloneClass )
    {
        task =
            CreateTaskByRegisteredFunction(
                cloneClass );
    }
    else
    {
        // Use the default task implementation.
        task = new SimpleTask();
    }

    if ( task == nullptr )
    {
        return nullptr;
    }

    task->taskId = operationId;
    task->taskName =
        MessageMap::GetMsgName( operationId );

    SolverState::solverType = solverType;

    SetTaskAction( task );

    return task;
}

// ============================================================
// Resource preparation
// ============================================================

void ConfigureTaskFile( Task * task, int operationId, int solverType )
{
    HXClone * cloneClass =
        ONEFLOW::GetClass(
            operationId,
            solverType,
            FILE_FUNC );

    if ( cloneClass )
    {
        ONEFLOW::ConfigureTaskFile(
            task,
            cloneClass->data );
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

void SingleSolverSingleGridTask( int operationId )
{
    // Runtime form: plan + execute by id (no string lookup here).
    GenerateCmdList( operationId );
    CMD::ExecuteCmd();
}

void SingleSolverSingleGridTask( const std::string & taskName )
{
    // Source form: resolve name once, then use id path.
    const int operationId = MessageMap::GetMsgId( taskName );
    SingleSolverSingleGridTask( operationId );
}


// ============================================================
// Multi-solver / multi-grid execution
// ============================================================

void MultiSolverMultiGridTask( int operationId )
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

            // Id path: no per-level MessageMap::GetMsgId.
            ONEFLOW::SingleSolverSingleGridTask( operationId );
        }
    }
}

void MultiSolverMultiGridTask( const std::string & taskName )
{
    // Resolve once outside solver x grid loops.
    const int operationId = MessageMap::GetMsgId( taskName );
    MultiSolverMultiGridTask( operationId );
}



EndNameSpace