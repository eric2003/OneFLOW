/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class Task;

class Command
{
public:
    typedef HXVector< Task * > TList;
public:
    Command();
    virtual ~Command();
public:
    virtual void Execute() = 0;
public:
    TList * GetTaskList() { return tasks; };
    void AddTask( Task * task );
public:
    TList * tasks;
};

class NullCmd : public Command
{
public:
    NullCmd(){};
    ~NullCmd() override {};
public:
    void Execute() override {};
};

class SimpleCmd : public Command
{
public:
    SimpleCmd();
    ~SimpleCmd() override;
public:
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
    static void AddCmd( Command * cmd );
    static void RunCmd( Command * cmd );
    static void Clear();
    static void ExecuteCmd();
    static void ShowCmdInfo( Command * cmd, int iCmd );
public:
    static HXVector< Command * > * cmdList;
};

EndNameSpace
