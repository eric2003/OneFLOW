/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "UINsCom.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "SolverDef.h"
#include "Ctrl.h"
#include "TimeIntegral.h"

BeginNameSpace( ONEFLOW )

UINsField  uinsf;

UINsField::UINsField()
{
    ;
}

UINsField::~UINsField()
{
    ;
}

void UINsField::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    // INS mean-flow view - pin INC_NS_SOLVER (do not depend on ambient SolverState).
    const int ins = INC_NS_SOLVER;
    q  = FieldHome::GetSolutionField( grid, ins );
    q1  = FieldHome::GetStageSolutionField( grid, ins, 0 );
    q2  = FieldHome::GetStageSolutionField( grid, ins, 1 );
    dq  = FieldHome::GetUpdateDeltaField( grid, ins );
    limiter  = GetFieldPointer< MRField > ( grid, "limiter" );
    gama  = GetFieldPointer< MRField > ( grid, "gama" );
    dqdx  = FieldHome::GetGradField( grid, ins, 1 );
    dqdy  = FieldHome::GetGradField( grid, ins, 2 );
    dqdz  = FieldHome::GetGradField( grid, ins, 3 );

    dtdx  = GetFieldPointer< MRField > ( grid, "dtdx" );
    dtdy  = GetFieldPointer< MRField > ( grid, "dtdy" );
    dtdz  = GetFieldPointer< MRField > ( grid, "dtdz" );
    bc_q    = GetFieldPointer< MRField > ( grid, "bc_q" );
    bcdqdx  = FieldHome::GetBcGradField( grid, ins, 0 );
    bcdqdy  = FieldHome::GetBcGradField( grid, ins, 1 );
    bcdqdz  = FieldHome::GetBcGradField( grid, ins, 2 );

    visl  = GetFieldPointer< MRField > ( grid, "visl" );
    vist  = GetFieldPointer< MRField > ( grid, "vist" );

    tempr = GetFieldPointer< MRField >( grid, "tempr" );

    timestep = GetFieldPointer< MRField >( grid, "timestep" );

    invsr = ONEFLOW::GetFieldPointer< MRField >( grid, "invsr" );
    vissr = ONEFLOW::GetFieldPointer< MRField >( grid, "vissr" );
    impsr = ONEFLOW::GetFieldPointer< MRField >( grid, "impsr" );

    res  = FieldHome::GetResidualField( grid, ins );
    res1  = FieldHome::GetStageResidualField( grid, ins, 0 );
    res2  = FieldHome::GetStageResidualField( grid, ins, 1 );

    rhs = res;

    if ( ctrl.time_integral == MULTI_STAGE )
    {
        // Alias: update delta shares residual storage for multi-stage schemes.
        dq = res;
    }
}

EndNameSpace
