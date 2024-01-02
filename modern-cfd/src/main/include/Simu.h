/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include <json/json.h>

//Types of problem solving
enum class TaskLineEnum
{
    UNDEFINED = 0,
    SOLVE_FIELD = 1,
    GRID_GEN = 2
};

class Solver_t
{
public:
    Solver_t() {};
    ~Solver_t() {};
public:
    virtual void Init( Json::Value & root ) {};
    virtual void Run() {};
};

class CfdPara;
class FieldSolver_t : public Solver_t
{
public:
    FieldSolver_t();
    ~FieldSolver_t();
public:
    void Init( Json::Value & root ) override;
    void Run() override;
private:
    CfdPara * cfd_para;
};

class GridSolver_t : public Solver_t
{
public:
    GridSolver_t() {};
    ~GridSolver_t() {};
public:
    void Init( Json::Value & root ) override;
    void Run() override;
private:
    int gridobj;
    std::string gridName;
};


class Simu
{
public:
    Simu( int argc, char ** argv );
    ~Simu();
public:
    void Init( int argc, char ** argv );
    void Run();
    void ReadControlParameter();
    void Process( Json::Value & root );
public:
    TaskLineEnum GetTaskLine();
    Solver_t * solver_t;
};

