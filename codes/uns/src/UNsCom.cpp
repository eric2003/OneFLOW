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

#include "UCom.h"
#include "UNsCom.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"
#include "Ctrl.h"
#include "TimeIntegral.h"

BeginNameSpace( ONEFLOW )

UNsField unsf;

UNsField::UNsField()
{
    ;
}

UNsField::~UNsField()
{
    ;
}

void UNsField::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    q  = GetFieldPointer< MRField > ( grid, "q" );
    q1  = GetFieldPointer< MRField > ( grid, "q1" );
    q2  = GetFieldPointer< MRField > ( grid, "q2" );
    dq  = GetFieldPointer< MRField > ( grid, "dq" );
    limiter  = GetFieldPointer< MRField > ( grid, "limiter" );
    gama  = GetFieldPointer< MRField > ( grid, "gama" );
    dqdx  = GetFieldPointer< MRField > ( grid, "dqdx" );
    dqdy  = GetFieldPointer< MRField > ( grid, "dqdy" );
    dqdz  = GetFieldPointer< MRField > ( grid, "dqdz" );

    dtdx  = GetFieldPointer< MRField > ( grid, "dtdx" );
    dtdy  = GetFieldPointer< MRField > ( grid, "dtdy" );
    dtdz  = GetFieldPointer< MRField > ( grid, "dtdz" );
    bc_q    = GetFieldPointer< MRField > ( grid, "bc_q" );
    bcdqdx  = GetFieldPointer< MRField > ( grid, "bcdqdx" );
    bcdqdy  = GetFieldPointer< MRField > ( grid, "bcdqdy" );
    bcdqdz  = GetFieldPointer< MRField > ( grid, "bcdqdz" );

    visl  = GetFieldPointer< MRField > ( grid, "visl" );
    vist  = GetFieldPointer< MRField > ( grid, "vist" );

    tempr = GetFieldPointer< MRField >( grid, "tempr" );

    timestep = GetFieldPointer< MRField >( grid, "timestep" );

    invsr = ONEFLOW::GetFieldPointer< MRField >( grid, "invsr" );
    vissr = ONEFLOW::GetFieldPointer< MRField >( grid, "vissr" );
    impsr = ONEFLOW::GetFieldPointer< MRField >( grid, "impsr" );

    res  = GetFieldPointer< MRField > ( grid, "res" );
    res1  = GetFieldPointer< MRField > ( grid, "res1" );
    res2  = GetFieldPointer< MRField > ( grid, "res2" );

    rhs = res;

    if ( ctrl.time_integral == MULTI_STAGE )
    {
        dq = res;
    }
}

EndNameSpace