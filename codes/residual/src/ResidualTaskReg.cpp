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
#include "ResidualTaskReg.h"
#include "ResidualTask.h"
#include "SolverState.h"
#include "UResidual.h"

BeginNameSpace( ONEFLOW )

void RegisterRedisualTask()
{
    REGISTER_DATA_CLASS( CreateResidualTask );
    REGISTER_DATA_CLASS( DumpResidual );
}

void CreateResidualTask( StringField & data )
{
    cout << "CreateResidualTask( StringField & data )\n";
    ResidualTask * task = new ResidualTask();
    cout << "CreateResidualTask 1\n";
    TaskState::task = task;
    cout << "CreateResidualTask 2\n";
}

void DumpResidual( StringField & data )
{
    cout << " DumpResidual ( StringField & data ) \n";
    int sTid = SolverState::tid;

    cout << " sTid = " << sTid << "\n";

    Residual * residual = new UResidual();
    cout << " DumpResidual 1 \n";
    residual->Dump( sTid );
    cout << " DumpResidual 2 \n";
    delete residual;
    cout << " DumpResidual 3 \n";
}

EndNameSpace