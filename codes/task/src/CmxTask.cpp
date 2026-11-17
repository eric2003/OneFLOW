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
#include "StageProfiler.h"
#include "FileMap.h"
#include <memory>
#include <utility>
#include <map>
#include <cstdint>
#include <chrono>

BeginNameSpace( ONEFLOW )

namespace {

// Times one plan+execute dispatch.  When the stage breakdown file is unset the
// object is inert, so the production path keeps paying only the Enabled()
// environment lookup that the coarse-stage timers already perform.
struct TaskDispatchTimer
{
    bool active = false;
    std::chrono::steady_clock::time_point start;

    TaskDispatchTimer() : active( StageProfiler::Enabled() )
    {
        if ( active ) start = std::chrono::steady_clock::now();
    }

    void Finish( int operationId, const char * taskName ) const
    {
        if ( ! active ) return;
        const auto end = std::chrono::steady_clock::now();
        const std::chrono::duration< double, std::milli > elapsed = end - start;
        StageProfiler::RecordTask( operationId, taskName, elapsed.count() );
    }
};
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

std::unique_ptr<Task> CreateTaskByRegisteredFunction(
    HXClone * cloneClass )
{
    if ( cloneClass == nullptr )
    {
        return nullptr;
    }

    // TASK_FUNC callbacks place their construction result here.
    TaskState::createdTask.reset();

    cloneClass->Solve();

    // Take ownership; leave the slot empty.
    return std::move( TaskState::createdTask );
}

const char * GetFunctionTypeName(
    int funcType )
{
    switch ( funcType )
    {
    case COMM_FUNC:
        return "COMM_FUNC";

    case RECV_FUNC:
        return "RECV_FUNC";

    case MESG_FUNC:
        return "MESG_FUNC";

    case TASK_FUNC:
        return "TASK_FUNC";

    case FILE_FUNC:
        return "FILE_FUNC";

    default:
        return "UNKNOWN";
    }
}

bool HasRegisteredFunction(
    int operationId,
    int solverType,
    int funcType )
{
    HXRegister * hxRegister =
        RegisterFactory::GetRegister(
            solverType,
            funcType );

    if ( hxRegister == nullptr )
    {
        return false;
    }

    const std::string & operationName =
        MessageMap::GetMsgName( operationId );

    return hxRegister->GetClass(
        operationName ) != nullptr;
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

    // The register/class lookup is epoch-cached by GetClass().  The old path
    // repeated the MessageMap name lookup and HXRegister::GetClass string map
    // on every RK task (thousands of calls per run); reuse the cache.
    HXClone * cloneClass = GetClass( operationId, solverType, MESG_FUNC );

    if ( cloneClass )
    {
        // Expand the operation into runtime commands.
        cloneClass->Solve();
    }
    else
    {
        // Treat the operation as a single runtime command.
        ONEFLOW::AddCmdToList( MessageMap::GetMsgName( operationId ) );
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
    // Ownership pipeline (unique_ptr only):
    //   CreateTask -> ConfigureTaskFile(task.get()) -> SimpleCmd::AddTask(move)
    //   -> CMD::AddCmd(move). No raw new/delete on the hot path.
    auto task =
        ONEFLOW::CreateTask(
            operationId,
            solverType );

    if ( task == nullptr )
    {
        return;
    }

    // Prepare files/resources required by the operation.
    ONEFLOW::ConfigureTaskFile(
        task.get(),
        operationId,
        solverType );

    // Build the command with RAII ownership.
    auto cmd = std::make_unique<SimpleCmd>();

    // Transfer Task ownership to the Command.
    cmd->AddTask( std::move( task ) );

    // Transfer Command ownership to CMD.
    CMD::AddCmd( std::move( cmd ) );
}


// ============================================================
// Task construction
// ============================================================

std::unique_ptr<Task> CreateTask( int operationId, int solverType )
{
    std::unique_ptr<Task> task;

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
        task = std::make_unique<SimpleTask>();
    }

    if ( task == nullptr )
    {
        return nullptr;
    }

    task->taskId = operationId;
    task->taskName =
        MessageMap::GetMsgName( operationId );

    SolverState::solverType = solverType;

    SetTaskAction( task.get() );

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

 /*   if ( funcType == COMM_FUNC &&
        ( TaskState::task->taskName == "INIT_FIRST" ||
        TaskState::task->taskName == "INIT_RESTART" ||
        TaskState::task->taskName == "READ_RESTART" ||
        TaskState::task->taskName == "INIT_FINAL" ) 
        )
    {
        DumpTaskDispatch( std::cout );
    }*/


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

void DumpTaskDispatch(
    std::ostream & output )
{
    if ( TaskState::task == nullptr )
    {
        output
            << "Task Dispatch: <no current task>\n";
        return;
    }

    const int operationId =
        TaskState::task->taskId;

    const int solverType =
        SolverState::solverType;

    output
        << "========== Task Dispatch Environment ==========\n\n";

    output
        << "Task:\n"
        << "  operation : "
        << TaskState::task->taskName
        << '\n'
        << "  taskId    : "
        << operationId
        << '\n'
        << "  solverType: "
        << solverType
        << "\n\n";

    const int functionTypes[] =
    {
        MESG_FUNC,
        TASK_FUNC,
        FILE_FUNC,
        COMM_FUNC,
        RECV_FUNC
    };

    for ( int i = 0; i < 5; ++ i )
    {
        const int funcType =
            functionTypes[ i ];

        output
            << "  "
            << GetFunctionTypeName( funcType )
            << " : "
            << (
                HasRegisteredFunction(
                    operationId,
                    solverType,
                    funcType )
                ? "registered"
                : "none"
                )
            << '\n';
    }

    output
        << "\n===============================================\n";
}

// ============================================================
// Operation execution entry
// ============================================================

void SingleSolverSingleGridTask( int operationId )
{
    TaskDispatchTimer dispatchTimer;
    const std::string & taskName = MessageMap::GetMsgName( operationId );

    // Runtime form: plan + execute by id (no string lookup here).
    {
        ScopedStageTimer planTimer( "task_plan" );
        GenerateCmdList( operationId );
    }
    const char * category = StageProfiler::CategoryForTask( taskName );
    if ( category[ 0 ] != '\0' )
    {
        ScopedStageTimer timer( category );
        CMD::ExecuteCmd();
    }
    else
    {
        CMD::ExecuteCmd();
    }
    dispatchTimer.Finish( operationId, taskName.c_str() );
}

void SingleSolverSingleGridTask(
    int operationId,
    SimuContext & context )
{
    TaskDispatchTimer dispatchTimer;
    const std::string & taskName = MessageMap::GetMsgName( operationId );

    {
        ScopedStageTimer planTimer( "task_plan" );
        GenerateCmdList( operationId );
    }
    const char * category =
        StageProfiler::CategoryForTask( taskName );
    if ( category[ 0 ] != '\0' )
    {
        ScopedStageTimer timer( category );
        CMD::ExecuteCmd( & context );
    }
    else
    {
        CMD::ExecuteCmd( & context );
    }
    dispatchTimer.Finish( operationId, taskName.c_str() );
}

void SingleSolverSingleGridTask( const std::string & taskName )
{
    // Source form: resolve name once, then use id path.
    const int operationId = MessageMap::GetMsgId( taskName );
    SingleSolverSingleGridTask( operationId );
}

void SingleSolverSingleGridTask(
    const std::string & taskName,
    SimuContext & context )
{
    const int operationId = MessageMap::GetMsgId( taskName );
    SingleSolverSingleGridTask( operationId, context );
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

void MultiSolverMultiGridTask(
    int operationId,
    SimuContext & context )
{
    for ( int solverIndex = 0;
        solverIndex < SolverState::nSolver;
        ++ solverIndex )
    {
        SolverState::SetSolverTypeBySolverIndex( solverIndex );
        for ( int gl = 0; gl < GridState::nGrids; ++ gl )
        {
            GridState::SetGridLevel( gl );
            SingleSolverSingleGridTask( operationId, context );
        }
    }
}

void MultiSolverMultiGridTask( const std::string & taskName )
{
    // Resolve once outside solver x grid loops.
    const int operationId = MessageMap::GetMsgId( taskName );
    MultiSolverMultiGridTask( operationId );
}

void MultiSolverMultiGridTask(
    const std::string & taskName,
    SimuContext & context )
{
    const int operationId = MessageMap::GetMsgId( taskName );
    MultiSolverMultiGridTask( operationId, context );
}


EndNameSpace
