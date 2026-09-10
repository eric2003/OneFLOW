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

BeginNameSpace( ONEFLOW )

class Task;

class Command
{
public:
    using TList = HXVector< Task * >;
    using TaskOwnerList = HXVector< std::unique_ptr< Task > >;

public:
    Command();
    virtual ~Command();

    Command( const Command & ) = delete;
    Command & operator=( const Command & ) = delete;

    // Moving is disabled because 'tasks' points to taskList_.
    Command( Command && ) = delete;
    Command & operator=( Command && ) = delete;

public:
    virtual void Execute() = 0;

public:
    TList * GetTaskList() { return tasks; }
    const TList * GetTaskList() const { return tasks; }

    // Takes ownership of the raw Task pointer.
    // This overload preserves compatibility with existing code.
    void AddTask( Task * task );

    // Preferred RAII interface for new code.
    void AddTask( std::unique_ptr< Task > task );

public:
    // Non-owning compatibility view.
    TList * tasks;

private:
    TList taskList_;
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
    SimpleCmd();
    ~SimpleCmd() override;

    void Execute() override;
};

class CMD
{
public:
    CMD();
    ~CMD();

public:
    static void Init();
    static void Free();

public:
    // Takes ownership of the raw Command pointer.
    // This overload preserves the existing API.
    static void AddCmd( Command * cmd );

    // Preferred RAII interface for new code.
    static void AddCmd( std::unique_ptr< Command > cmd );

    // Execute without transferring ownership.
    static void RunCmd( Command * cmd );

    // Destroy all queued commands.
    static void Clear();

    // Execute the current command batch and then destroy it.
    static void ExecuteCmd();

    static void ShowCmdInfo( Command * cmd, int iCmd );

public:
    // Non-owning compatibility view of the owned command queue.
    static HXVector< Command * > * cmdList;

private:
    using CommandOwnerList =
        HXVector< std::unique_ptr< Command > >;

    // Actual owner of all commands currently in the queue.
    static CommandOwnerList * commandOwners;
};

EndNameSpace