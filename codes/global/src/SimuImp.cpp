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
#include "GridFactory.h"
#include "ParaFile.h"
#include "Parallel.h"
#include "Test.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void RunSimu()
{
    simu_state.Init();

	//设置模拟参数
	if (simu_state.simutask != TaskEnum::FUN_TEST)
	{
		ConstructSystemMap();
	}

	//根据不同的simutask值，执行不同的求解流程
    if ( simu_state.simutask == TaskEnum::SOLVE_FIELD )
    {
        FieldSimu();
    }
    else if ( simu_state.simutask == TaskEnum::CREATE_GRID )
    {
        GenerateGrid();
    }
    else if ( simu_state.simutask == TaskEnum::CREATE_WALL_DIST )
    {
        WalldistSimu();
    }
    else if ( simu_state.simutask == TaskEnum::FUN_TEST )
    {
        FunTest();
    }
}

void InitSimu()
{
	cout << "OneFlow running\n";
	ONEFLOW::SetUpParallelEnvironment();
	ONEFLOW::ReadControlInformation();
}

EndNameSpace