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

#include "UTurbGrad.h"
#include "UCom.h"
#include "TurbCom.h"
#include "DataBase.h"
#include "FieldImp.h"
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

    name  = "turbq";
    namex = "turbdqdx";
    namey = "turbdqdy";
    namez = "turbdqdz";

    q    = GetFieldPointer< MRField > ( grid, name  );
    dqdx = GetFieldPointer< MRField > ( grid, namex );
    dqdy = GetFieldPointer< MRField > ( grid, namey );
    dqdz = GetFieldPointer< MRField > ( grid, namez );

    bdqdx = GetFieldPointer< MRField > ( grid, "turbbcdqdx" );
    bdqdy = GetFieldPointer< MRField > ( grid, "turbbcdqdy" );
    bdqdz = GetFieldPointer< MRField > ( grid, "turbbcdqdz" );

    this->nEqu = turbcom.nEqu;

    this->istore = 1;
}

EndNameSpace