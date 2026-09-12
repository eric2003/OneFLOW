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
#include "Command.h"
#include "Task.h"
#include "TaskState.h"

#include <iostream>
#include <memory>
#include <string>
#include <utility>

BeginNameSpace( ONEFLOW )

/*
* Command implementation
*
* The actual Task ownership is maintained by ownedTasks_.
* taskList_ only provides a non-owning compatibility view.
*/

Command::Command() = default;
Command::~Command() = default;

void Command::AddTask( Task * task )
{
    if ( task == nullptr )
    {
        return;
    }

    // Preserve the legacy raw-pointer interface while transferring
    // ownership to the RAII implementation.
    AddTask( std::unique_ptr< Task >( task ) );
}

void Command::AddTask( std::unique_ptr< Task > task )
{
    if ( ! task )
    {
        return;
    }

    Task * rawTask = task.get();

    // Transfer ownership to the Command.
    ownedTasks_.push_back( std::move( task ) );

    // Keep the execution-order view non-owning.
    taskList_.push_back( rawTask );
}

SimpleCmd::SimpleCmd() = default;

SimpleCmd::~SimpleCmd() = default;

void SimpleCmd::Execute()
{
    const TList & taskList = * GetTaskList();

    for ( HXSize_t iTask = 0; iTask < taskList.size(); ++ iTask )
    {
        Task * task = taskList[ iTask ];

        if ( task == nullptr )
        {
            continue;
        }

        // Set the currently executing Task.
        TaskState::task = task;

        task->Run();
    }
}


/*
* CMD implementation
*
* commandOwners is the actual owner of every Command in the queue.
* cmdList is only a non-owning compatibility view.
*/

HXVector< Command * > * CMD::cmdList_ = nullptr;
CMD::CommandOwnerList * CMD::commandOwners = nullptr;

CMD::CMD()
{
}

CMD::~CMD()
{
}

void CMD::Init()
{
    if ( CMD::cmdList_ != 0 )
    {
        return;
    }

    CMD::cmdList_ = new HXVector< Command * >;
    CMD::commandOwners = new CommandOwnerList;
}

void CMD::Free()
{
    /*
    * Clear() releases every owned Command first.
    *
    * This is important because simply deleting cmdList_ would only
    * destroy the non-owning pointer container.
    */
    CMD::Clear();

    delete CMD::commandOwners;
    CMD::commandOwners = nullptr;

    delete CMD::cmdList_;
    CMD::cmdList_ = nullptr;

    TaskState::task = nullptr;
}

//void CMD::AddCmd( Command * cmd )
//{
//    if ( cmd == nullptr )
//    {
//        return;
//    }
//
//    CMD::Init();
//
//    /*
//    * The raw pointer API is retained for compatibility.
//    * Ownership is transferred immediately to commandOwners.
//    */
//    CMD::cmdList_->push_back( cmd );
//    CMD::commandOwners->emplace_back( cmd );
//}

void CMD::AddCmd( Command * cmd )
{
    if ( cmd == nullptr )
    {
        return;
    }

    // Preserve the legacy raw-pointer interface while transferring
    // ownership to the RAII implementation.
    CMD::AddCmd( std::unique_ptr< Command >( cmd ) );
}

void CMD::AddCmd( std::unique_ptr< Command > cmd )
{
    //if ( ! cmd )
    //{
    //    return;
    //}

    if ( cmd == nullptr )
    {
        return;
    }

    CMD::Init();

    Command * rawCmd = cmd.get();

    /*
    * Transfer ownership into the queue.
    */
    CMD::commandOwners->push_back( std::move( cmd ) );

    /*
    * Keep the non-owning execution view synchronized.
    */
    CMD::cmdList_->push_back( rawCmd );
}

const HXVector< Command * > * CMD::GetCmdList()
{
    return CMD::cmdList_;
}

void CMD::RunCmd( Command * cmd )
{
    if ( cmd == nullptr )
    {
        return;
    }

    cmd->Execute();
}

void CMD::Clear()
{
    if ( CMD::cmdList_ == nullptr )
    {
        return;
    }

    /*
    * Destroy Commands first.
    *
    * Command destruction also destroys all Tasks owned by the Command.
    */
    if ( CMD::commandOwners != nullptr )
    {
        CMD::commandOwners->clear();
    }

    /*
    * cmdList_ contains only non-owning pointers.
    * Clear it after destroying the actual objects so that no stale
    * pointers remain visible through the read-only compatibility view.
    */
    CMD::cmdList_->clear();

    TaskState::task = nullptr;
}

void CMD::ExecuteCmd()
{
    if ( CMD::cmdList_ == nullptr )
    {
        return;
    }

    /*
    * Keep the original batch semantics:
    *
    * - Only Commands that existed when ExecuteCmd() started are
    *   executed in this dispatch cycle.
    * - Commands are executed in FIFO order.
    *
    * Commands added while a Command is executing are not executed
    * in the current batch. Clear() subsequently releases them.
    */
    const HXSize_t nCmd = CMD::cmdList_->size();

    for ( HXSize_t iCmd = 0; iCmd < nCmd; ++ iCmd )
    {
        Command * cmd = ( * CMD::cmdList_ )[ iCmd ];

        if ( cmd == nullptr )
        {
            continue;
        }

        // CMD::ShowCmdInfo( cmd, static_cast< int >( iCmd ) );

        cmd->Execute();
    }

    /*
    * No individual delete is required here.
    *
    * commandOwners owns all Commands and Clear() releases them
    * through unique_ptr.
    */
    CMD::Clear();
}

void CMD::ShowCmdInfo( Command * cmd, int iCmd )
{
    if ( cmd == nullptr )
    {
        return;
    }

    const Command::TList & tasks = * cmd->GetTaskList();

    for ( HXSize_t i = 0; i < tasks.size(); ++ i )
    {
        Task * task = tasks[ i ];

        if ( task == nullptr )
        {
            continue;
        }

        int iTaskGlobal =
            iCmd + static_cast< int >( i );

        std::cout
            << " iTaskGlobal = " << iTaskGlobal
            << " iTaskLocal = " << i
            << " TaskCode = " << task->taskId
            << " Task Name = " << task->taskName
            << std::endl;
    }
}

EndNameSpace