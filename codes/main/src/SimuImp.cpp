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
#include "SimuImp.h"
#include "SimuDef.h"
#include "SimpleSimu.h"
#include "System.h"
#include "FieldSimu.h"
#include "MultiBlock.h"
#include "Prj.h"
#include "ParaFile.h"
#include "Parallel.h"
#include "GridFactory.h"
#include "Test.h"
#include "Theory.h"
#include "PostProcess.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

SimuImp::SimuImp( std::vector<std::string> &args )
{
    this->args = args;
    Prj::ProcessCmdLineArgs( args );
}

SimuImp::~SimuImp()
{
}

void SimuImp::Run()
{
    this->PreProcess();
    this->MainProcess();
    this->PostProcess();
}

void SimuImp::PreProcess()
{
    InitSimu();
}

void SimuImp::MainProcess()
{
    RunSimu();
}

void SimuImp::PostProcess()
{
    HXFinalize();
}

void SimuImp::RunSimu()
{
    //Set the type of operation that ONEFLOW needs to perform
    simu_state.Init();

    //Call different solving modules according to the task type
    const TaskEnum task = simu_state.Task();

    if ( task == TaskEnum::SOLVE_FIELD ||
         task == TaskEnum::CREATE_GRID ||
         task == TaskEnum::CREATE_WALL_DIST
       )
    {
        ConstructSystemMap();
    }

    //According to different simutask values, different solving processes are executed
    switch ( task )
    {
        case TaskEnum::SOLVE_FIELD:
            FieldSimu();
            break;
        case TaskEnum::CREATE_GRID:
            GenerateGrid();
            break;
        case TaskEnum::CREATE_WALL_DIST:
            WalldistSimu();
            break;
        case TaskEnum::FUNCTION_TEST:
            FunctionTest();
            break;
        case TaskEnum::SOLVE_THEORY:
            TheorySimu();
            break;
        case TaskEnum::TOY_MODEL:
            ToyModelSimu();
            break;
        case TaskEnum::POST_TASK:
            PostSimu();
            break;
        default:
        {
            std::cerr << "unknown simutask value!!" << std::endl;
            exit(EXIT_FAILURE);
        }
        break;
    }
}

void SimuImp::InitSimu()
{
    std::cout << " OneFLOW is running\n";
    ONEFLOW::SetUpParallelEnvironment();
    ONEFLOW::ReadControlInfo();
}


EndNameSpace
