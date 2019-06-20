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
#include "Simulation.h"
#include "SimuDef.h"
#include "System.h"
#include "Parallel.h"
#include "ParaFile.h"
#include "GridFactory.h"
#include "FieldSimu.h"
#include "MultiBlock.h"
#include "Test.h"
#include "Prj.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Simulation::Simulation( int argc, char ** argv )
{
	args.resize(argc);

	cout << " args.size() = " << args.size() << "\n";

    for ( int i = 0; i < argc; ++ i )
    {   
		args[i] = argv[i];
        cout<<"arguments["<<i<<"] is: " << args[ i ] << endl;
    }

	if ( args.size() > 1 )
	{
		string prjName = args[1];
		PrjStatus::SetPrjBaseDir( prjName );
	}
}

Simulation::~Simulation()
{
	if (!args.empty()) 
	{
		args.clear();
		args.shrink_to_fit();
	};
}

void Simulation::Run()
{
	this->PreProcess();
	this->MainProcess();
	this->PostProcess();
}

void Simulation::PreProcess()
{
    InitSimu();
}

void Simulation::MainProcess()
{
    RunSimu();
}

void Simulation::PostProcess()
{
}

void Simulation::RunSimu()
{
	//设置oneflow需要执行的操作类型
	simu_state.Init();

	//根据任务类型调用不同的求解模块
	const TaskEnum task = simu_state.Task();

	if (task != TaskEnum::FUN_TEST)
	{
		ConstructSystemMap();
	}

	//根据不同的simutask值，执行不同的求解流程
	switch (task)
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

void Simulation::InitSimu()
{
	cout << "OneFlow running\n";
	ONEFLOW::SetUpParallelEnvironment();
	ONEFLOW::ReadControlInformation();
}

EndNameSpace
