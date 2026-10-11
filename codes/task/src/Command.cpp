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

#include "DataBook.h"
#include "FileInfo.h"

#include <iostream>
#include <memory>
#include <string>
#include <utility>

BeginNameSpace( ONEFLOW )

namespace
{

    void DumpTaskInfo(
        std::ostream & output,
        Task * task,
        int iTask )
    {
        output
            << "  Task "
            << iTask
            << ":\n";

        if ( task == nullptr )
        {
            output
                << "    <null>\n";
            return;
        }

        output
            << "    id       : "
            << task->taskId
            << "\n"
            << "    name     : "
            << task->taskName
            << "\n";

        output
            << "    action:\n"
            << "      normal  : "
            << ( task->action != nullptr ? "set" : "null" )
            << "\n"
            << "      send    : "
            << ( task->sendAction != nullptr ? "set" : "null" )
            << "\n"
            << "      receive : "
            << ( task->recvAction != nullptr ? "set" : "null" )
            << "\n";

        if ( task->dataBook != nullptr )
        {
            output
                << "    DataBook:\n"
                << "      size      : "
                << task->dataBook->size()
                << "\n"
                << "      pageCount : "
                << task->dataBook->GetPageCount()
                << "\n";
        }
        else
        {
            output
                << "    DataBook: <null>\n";
        }

        if ( task->fileInfo != nullptr )
        {
            output
                << "    FileInfo:\n"
                << "      fileName : "
                << task->fileInfo->fileName
                << "\n";
        }
        else
        {
            output
                << "    FileInfo: <null>\n";
        }
    }

}

// Command: single ownership vector ownedTasks_ (order == AddTask order).

void Command::AddTask( Task * task )
{
    if ( task == nullptr )
    {
        return;
    }

    AddTask( std::unique_ptr< Task >( task ) );
}

void Command::AddTask( std::unique_ptr< Task > task )
{
    if ( ! task )
    {
        return;
    }

    ownedTasks_.push_back( std::move( task ) );
}

void SimpleCmd::Execute()
{
    const HXSize_t nTask = GetTaskCount();

    for ( HXSize_t iTask = 0; iTask < nTask; ++ iTask )
    {
        Task * task = GetTask( iTask );

        if ( task == nullptr )
        {
            continue;
        }

        TaskState::task = task;
        task->Run();
    }
}


// CMD: single ownership vector commandOwners (order == AddCmd order).

std::unique_ptr< CMD::CommandOwnerList > CMD::commandOwners;

void CMD::Init()
{
    if ( CMD::commandOwners != nullptr )
    {
        return;
    }

    CMD::commandOwners = std::make_unique< CommandOwnerList >();
}

void CMD::Free()
{
    CMD::Clear();
    CMD::commandOwners.reset();
    TaskState::task = nullptr;
}

void CMD::AddCmd( Command * cmd )
{
    if ( cmd == nullptr )
    {
        return;
    }

    CMD::AddCmd( std::unique_ptr< Command >( cmd ) );
}

void CMD::AddCmd( std::unique_ptr< Command > cmd )
{
    if ( cmd == nullptr )
    {
        return;
    }

    CMD::Init();
    CMD::commandOwners->push_back( std::move( cmd ) );
}

HXSize_t CMD::GetCmdCount()
{
    if ( CMD::commandOwners == nullptr )
    {
        return 0;
    }
    return CMD::commandOwners->size();
}

Command * CMD::GetCmd( HXSize_t index )
{
    if ( CMD::commandOwners == nullptr )
    {
        return nullptr;
    }
    return ( * CMD::commandOwners )[ index ].get();
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
    if ( CMD::commandOwners == nullptr )
    {
        return;
    }

    // unique_ptr elements destroy Commands (and their Tasks) on clear.
    CMD::commandOwners->clear();
    TaskState::task = nullptr;
}

void CMD::ExecuteCmd()
{
    if ( CMD::commandOwners == nullptr )
    {
        return;
    }

    // Batch semantics: only commands present at entry are executed.
    // Commands enqueued during Execute run in a later batch after Clear.
    const HXSize_t nCmd = CMD::commandOwners->size();

    for ( HXSize_t iCmd = 0; iCmd < nCmd; ++ iCmd )
    {
        Command * cmd = ( * CMD::commandOwners )[ iCmd ].get();

        if ( cmd == nullptr )
        {
            continue;
        }

        cmd->Execute();
    }

    CMD::Clear();
}

void CMD::ShowCmdInfo( Command * cmd, int iCmd )
{
    if ( cmd == nullptr )
    {
        return;
    }

    const HXSize_t nTask = cmd->GetTaskCount();

    for ( HXSize_t i = 0; i < nTask; ++ i )
    {
        Task * task = cmd->GetTask( i );

        if ( task == nullptr )
        {
            continue;
        }

        const int iTaskGlobal = iCmd + static_cast< int >( i );

        std::cout
            << " iTaskGlobal = " << iTaskGlobal
            << " iTaskLocal = " << i
            << " TaskCode = " << task->taskId
            << " Task Name = " << task->taskName
            << std::endl;
    }
}

void CMD::DumpCommandQueue( std::ostream & output )
{
    output
        << "========== Command Environment ==========\n";

    if ( CMD::commandOwners == nullptr ||
        CMD::commandOwners->empty() )
    {
        output
            << "  <empty>\n"
            << "==========================================\n";
        return;
    }

    for ( HXSize_t iCmd = 0;
        iCmd < CMD::commandOwners->size();
        ++ iCmd )
    {
        Command * cmd =
            ( * CMD::commandOwners )[ iCmd ].get();

        output
            << "\n[Command "
            << iCmd
            << "]\n";

        if ( cmd == nullptr )
        {
            output
                << "  <null>\n";
            continue;
        }

        const HXSize_t nTask = cmd->GetTaskCount();

        output
            << "  Tasks: "
            << nTask
            << "\n";

        for ( HXSize_t iTask = 0; iTask < nTask; ++ iTask )
        {
            DumpTaskInfo(
                output,
                cmd->GetTask( iTask ),
                static_cast< int >( iTask ) );
        }
    }

    output
        << "\n==========================================\n";
}


EndNameSpace
