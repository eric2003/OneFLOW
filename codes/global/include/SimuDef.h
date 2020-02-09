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
#pragma once
#include <map>
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

//求解问题的类型
enum class TaskEnum
{
    SOLVE_FIELD = 0,
    CREATE_GRID = 1,
    CREATE_WALL_DIST = 2,
    PARTITION_GRID = 3,
    FUN_TEST = 4
};
const map<string, TaskEnum> TaskFilter = 
{
    {"Solve",TaskEnum::SOLVE_FIELD},
    {"Grid",TaskEnum::CREATE_GRID},
    {"WallDist",TaskEnum::CREATE_WALL_DIST},
    {"Partition",TaskEnum::PARTITION_GRID},
    {"FunTest",TaskEnum::FUN_TEST}
};


//管理oneflow求解时执行的任务类型
class SimuState
{
public:
    SimuState();
    virtual ~SimuState();
public:
    //根据数据库的参数，设置simutask对应的值
    void Init();
    //返回要执行的任务类型
    const TaskEnum Task() const;

private:
    TaskEnum simutask;
};

extern SimuState simu_state;

EndNameSpace