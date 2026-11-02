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

#pragma once

#include "HXDefine.h"
#include <memory>
#include <ostream>

BeginNameSpace( ONEFLOW )

class Task;
class SimuContext;

// Command owns its Tasks via unique_ptr. Execution order is the order of
// AddTask() calls (index order in the owner vector).
class Command
{
public:
    using TaskOwnerList = HXVector< std::unique_ptr< Task > >;

public:
    Command() = default;
    virtual ~Command() = default;

    Command( const Command & ) = delete;
    Command & operator=( const Command & ) = delete;
    Command( Command && ) = delete;
    Command & operator=( Command && ) = delete;

public:
    virtual void Execute() = 0;
    virtual void ExecuteWithContext( SimuContext * context )
    {
        (void) context;
        Execute();
    }

public:
    // Non-owning accessors (order == ownership vector order).
    HXSize_t GetTaskCount() const { return ownedTasks_.size(); }
    Task * GetTask( HXSize_t index ) const
    {
        return ownedTasks_[ index ].get();
    }

    // Takes ownership of the raw Task pointer (legacy compatibility).
    void AddTask( Task * task );

    // Preferred RAII interface for new code.
    void AddTask( std::unique_ptr< Task > task );

private:
    TaskOwnerList ownedTasks_;
};

class NullCmd : public Command
{
public:
    NullCmd() = default;
    ~NullCmd() override = default;

    void Execute() override {}
};

class SimpleCmd : public Command
{
public:
    SimpleCmd() = default;
    ~SimpleCmd() override = default;

    void Execute() override;
    void ExecuteWithContext( SimuContext * context ) override;
};

// Process-wide command queue. Owns every queued Command via unique_ptr.
// Execution order is the order of AddCmd() (index order in the owner vector).
// There is no parallel non-owning pointer list.
class CMD
{
public:
    CMD() = default;
    ~CMD() = default;

public:
    static void Init();
    static void Free();

public:
    // Takes ownership of the raw Command pointer (legacy compatibility).
    static void AddCmd( Command * cmd );

    // Preferred RAII interface for new code.
    static void AddCmd( std::unique_ptr< Command > cmd );

    // Execute without transferring ownership.
    static void RunCmd( Command * cmd, SimuContext * context = nullptr );

    // Destroy all queued commands.
    static void Clear();

    // Execute the current command batch and then destroy it.
    static void ExecuteCmd( SimuContext * context = nullptr );

    static void ShowCmdInfo( Command * cmd, int iCmd );

    // Dump the current command queue without changing execution state.
    static void DumpCommandQueue( std::ostream & output );

    // Non-owning accessors into the owned queue.
    // GetCmdCount() == 0 when uninitialized or empty after Clear/Free.
    static HXSize_t GetCmdCount();
    static Command * GetCmd( HXSize_t index );

private:
    using CommandOwnerList = HXVector< std::unique_ptr< Command > >;

    // Sole owner of all commands currently in the queue.
    static std::unique_ptr< CommandOwnerList > commandOwners;
};

EndNameSpace
