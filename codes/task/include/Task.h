/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include <ios>
using namespace std;

BeginNameSpace( ONEFLOW )

typedef void ( * TaskFunction )();

class DataBook;

class FileInfo
{
public:
    FileInfo();
    ~FileInfo();
public:
    string fileName;
    ios_base::openmode openMode;
};

class Task
{
public:
    Task();
    virtual ~Task();
public:
    virtual void Run(){};
public:
    int taskId;
    string taskName;
    TaskFunction action, sendAction, recvAction;
    DataBook * dataBook;
    FileInfo * fileInfo;
public:
};

class TaskState
{
public:
    TaskState();
    ~TaskState();
public:
    static Task * task;
};

class SimpleTask : public Task
{
public:
    SimpleTask() {};
    ~SimpleTask(){};
public:
    void Run();
protected:
};

EndNameSpace