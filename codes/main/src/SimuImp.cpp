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
#include "SimuImp.h"
#include "SimuDef.h"
#include "System.h"
#include "FieldSimu.h"
#include "MultiBlock.h"
#include "Prj.h"
#include "ParaFile.h"
#include "Parallel.h"
#include "GridFactory.h"
#include "Test.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

SimuImp::SimuImp( std::vector<std::string> &args )
{
    this->args = args;
    this->ProcessCmdLineArgs( args );
}

SimuImp::~SimuImp()
{
}

void SimuImp::ProcessCmdLineArgs( std::vector<std::string> &args )
{
    this->args = args;

    string prjName = args[ 2 ];
    PrjStatus::SetPrjBaseDir( prjName );
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
    //����oneflow��Ҫִ�еĲ�������
    simu_state.Init();

    //�����������͵��ò�ͬ�����ģ��
    const TaskEnum task = simu_state.Task();

    if ( task != TaskEnum::FUN_TEST )
    {
        ConstructSystemMap();
    }

    //���ݲ�ͬ��simutaskֵ��ִ�в�ͬ���������
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
        case TaskEnum::FUN_TEST:
            FunTest();
            break;
        default:
        {
            cerr << "unknown simutask value!!" << endl;
            exit(EXIT_FAILURE);
        }
        break;
    }
}

void SimuImp::InitSimu()
{
    cout << "OneFlow running\n";
    ONEFLOW::SetUpParallelEnvironment();
    ONEFLOW::ReadControlInfo();
}


EndNameSpace
