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

#include "UTurbCom.h"
#include "UnsGrid.h"
#include "CellMesh.h"
#include "Zone.h"
#include "Ctrl.h"
#include "TimeIntegral.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "ZoneSharedFields.h"
#include "MeanFlowView.h"
#include "SolverDef.h"

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
    // Own turb fields: pin TURB_SOLVER (Init may be called with ambiguous ambient state).
    const int turb = TURB_SOLVER;

    q   = FieldHome::GetSolutionField( grid, turb );
    q1  = FieldHome::GetStageSolutionField( grid, turb, 0 );
    q2  = FieldHome::GetStageSolutionField( grid, turb, 1 );
    dq  = FieldHome::GetUpdateDeltaField( grid, turb );

    dqdx = FieldHome::GetGradField( grid, turb, 1 );
    dqdy = FieldHome::GetGradField( grid, turb, 2 );
    dqdz = FieldHome::GetGradField( grid, turb, 3 );

    // Mean-flow coupling: explicit cross-solver view (NS only; non-owning).
    this->meanFlow.Bind( grid );

    {
        ZoneSharedFieldView shared;
        shared.Bind( grid );
        visl = shared.visl;
        vist = shared.vist;
        timestep = shared.timestep;
    }

    bld   = GetFieldPointer< MRField > ( grid, "bld" );

    res   = FieldHome::GetResidualField( grid, turb );
    impsr = GetFieldPointer< MRField > ( grid, "turbimpsr" );
    cross = GetFieldPointer< MRField > ( grid, "cross" );
    len_scale = GetFieldPointer< MRField > ( grid, "len_scale" );
    matrix_l = GetFieldPointer< MRField > ( grid, "turb_matrix_l" );
    matrix_r = GetFieldPointer< MRField > ( grid, "turb_matrix_r" );

    dist = & grid->GetCellMesh().dist;

    rhs = res;

    if ( ctrl.time_integral == MULTI_STAGE )
    {
        // Alias: update delta shares residual storage for multi-stage schemes.
        dq = res;
    }
}

EndNameSpace
