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

#include "UTurbGrad.h"
#include "UCom.h"
#include "TurbCom.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "FieldManager.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"


BeginNameSpace( ONEFLOW )

UTurbGrad uturb_grad;

UTurbGrad::UTurbGrad()
{
    ;
}

UTurbGrad::~UTurbGrad()
{
    ;
}

void UTurbGrad::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    // Current-solver solution/grad via SolverInfo::solutionGradNames.
    name  = FieldHome::GetGradName( 0 );
    namex = FieldHome::GetGradName( 1 );
    namey = FieldHome::GetGradName( 2 );
    namez = FieldHome::GetGradName( 3 );

    q    = FieldHome::GetSolutionField( grid );
    dqdx = FieldHome::GetGradField( grid, 1 );
    dqdy = FieldHome::GetGradField( grid, 2 );
    dqdz = FieldHome::GetGradField( grid, 3 );

    bdqdx = FieldHome::GetBcGradField( grid, 0 );
    bdqdy = FieldHome::GetBcGradField( grid, 1 );
    bdqdz = FieldHome::GetBcGradField( grid, 2 );

    this->nEqu = turbcom.nEqu;

    this->istore = 1;
}

EndNameSpace
