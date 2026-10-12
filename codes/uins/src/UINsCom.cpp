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

    q  = FieldHome::GetSolutionField( grid );  // INS -> "q"
    q1  = GetFieldPointer< MRField > ( grid, "q1" );
    q2  = GetFieldPointer< MRField > ( grid, "q2" );
    // dq: literal "dq". MULTI_STAGE rebinds dq = res (see below).
    dq  = GetFieldPointer< MRField > ( grid, "dq" );
    limiter  = GetFieldPointer< MRField > ( grid, "limiter" );
    gama  = GetFieldPointer< MRField > ( grid, "gama" );
    dqdx  = FieldHome::GetGradField( grid, 1 );
    dqdy  = FieldHome::GetGradField( grid, 2 );
    dqdz  = FieldHome::GetGradField( grid, 3 );

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

    res  = FieldHome::GetResidualField( grid );
    res1  = GetFieldPointer< MRField > ( grid, "res1" );
    res2  = GetFieldPointer< MRField > ( grid, "res2" );

    rhs = res;

    if ( ctrl.time_integral == MULTI_STAGE )
    {
        // Alias: update delta shares residual storage for multi-stage schemes.
        dq = res;
    }
}

EndNameSpace
