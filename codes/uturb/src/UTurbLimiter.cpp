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

#include "UTurbLimiter.h"
#include "TurbCom.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "DataBase.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

TurbLimField::TurbLimField()
{
    this->nEqu = turbcom.nEqu;
}

TurbLimField::~TurbLimField()
{
}

void TurbLimField::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    q       = GetFieldPointer< MRField > ( grid, "turbq" );
    dqdx    = GetFieldPointer< MRField > ( grid, "turbdqdx" );
    dqdy    = GetFieldPointer< MRField > ( grid, "turbdqdy" );
    dqdz    = GetFieldPointer< MRField > ( grid, "turbdqdz" );
    limiter = GetFieldPointer< MRField > ( grid, "turblimiter" );

    this->nEqu = q->GetNEqu();

    qf1 = std::make_unique<MRField>( this->nEqu, grid->nFaces );
    qf2 = std::make_unique<MRField>( this->nEqu, grid->nFaces );

    this->ckfun = & NoCheck;
}

TurbLimiter::TurbLimiter()
{
    limf = std::make_unique<TurbLimField>();
    limflag = ILMT_ZERO;
}

TurbLimiter::~TurbLimiter()
{
}

EndNameSpace
