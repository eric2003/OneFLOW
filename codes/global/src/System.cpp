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
#include "System.h"
#include "DimensionImp.h"
#include "ResidualTaskReg.h"
#include "AeroForceTaskReg.h"
#include "HeatFluxTaskReg.h"
#include "MultigridTaskReg.h"
#include "UnsteadyTaskReg.h"
#include "RestartTaskReg.h"
#include "ImplicitTaskReg.h"
#include "InterfaceTaskReg.h"
#include "VisualTaskReg.h"
#include "UpdateTaskReg.h"
#include "LhsTaskReg.h"
#include "FieldTaskReg.h"
#include "SolverTaskReg.h"
#include "FileMap.h"
#include "TaskImp.h"
#include "SolverDef.h"
#include "GridTask.h"
#include "NsSolverImp.h"
#include "TurbSolverImp.h"
#include "TaskRegister.h"

BeginNameSpace( ONEFLOW )

void ConstructSystemMap()
{
    ONEFLOW::SetDimension();
    TaskRegister::Run();
    //ONEFLOW::RegisterComTask();
    ONEFLOW::RegisterFileTask();
    ONEFLOW::RegisterRedisualTask();
    ONEFLOW::RegisterForceTask();
    ONEFLOW::RegisterHeatFluxTask();
    ONEFLOW::RegisterMultigridTask();
    ONEFLOW::RegisterUnsteadyTask();
    ONEFLOW::RegisterRestartTask();
    ONEFLOW::RegisterImplicitTask();
    ONEFLOW::RegisterInterfaceTask();
    ONEFLOW::RegisterVisualTask();
    ONEFLOW::RegisterUpdateTask();
    ONEFLOW::RegisterLhsTask();
    ONEFLOW::RegisterFieldTask();
    ONEFLOW::CreateSysMap();

    CreateMsgMap();

    HXVector< RegData * > regDataArray;
    GetRegData( regDataArray );

    RegisterSolverTask( regDataArray );
}

void GetRegData( HXVector< RegData * > & regDataArray )
{
    regDataArray.push_back( GetGridReg() );
    regDataArray.push_back( GetNsReg() );
    regDataArray.push_back( GetTurbReg() );
}

EndNameSpace