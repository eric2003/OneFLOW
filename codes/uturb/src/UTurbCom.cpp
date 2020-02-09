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

#include "UTurbCom.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "Zone.h"
#include "Ctrl.h"
#include "TimeIntegral.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )

UTurbField uturbf;

UTurbField::UTurbField()
{

}

UTurbField::~UTurbField()
{

}

void UTurbField::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    q   = GetFieldPointer< MRField > ( grid, "turbq" );
    q1  = GetFieldPointer< MRField > ( grid, "turbq1" );
    q2  = GetFieldPointer< MRField > ( grid, "turbq2" );
    dq  = GetFieldPointer< MRField > ( grid, "turbdq" );

    dqdx = GetFieldPointer< MRField > ( grid, "turbdqdx" );
    dqdy = GetFieldPointer< MRField > ( grid, "turbdqdy" );
    dqdz = GetFieldPointer< MRField > ( grid, "turbdqdz" );

    q_ns    = GetFieldPointer< MRField > ( grid, "q" );
    dqdx_ns = GetFieldPointer< MRField > ( grid, "dqdx" );
    dqdy_ns = GetFieldPointer< MRField > ( grid, "dqdy" );
    dqdz_ns = GetFieldPointer< MRField > ( grid, "dqdz" );

    visl  = GetFieldPointer< MRField > ( grid, "visl" );
    vist  = GetFieldPointer< MRField > ( grid, "vist" );

    bld   = GetFieldPointer< MRField > ( grid, "bld" );

    res   = GetFieldPointer< MRField > ( grid, "turbres" );
    impsr = GetFieldPointer< MRField > ( grid, "turbimpsr" );
    cross = GetFieldPointer< MRField > ( grid, "cross" );
    len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );
    matrix_l = GetFieldPointer< MRField > ( grid, "turb_matrix_l" );
    matrix_r = GetFieldPointer< MRField > ( grid, "turb_matrix_r" );
    timestep = GetFieldPointer< MRField > ( grid, "timestep" );

    dist = & grid->cellMesh->dist;

    rhs = res;

    if ( ctrl.time_integral == MULTI_STAGE )
    {
        dq = res;
    }
}

EndNameSpace