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

#include "Unsteady.h"
#include "UsdPara.h"
#include "FieldImp.h"
#include "FieldWrap.h"
#include "DataBase.h"
#include "Zone.h"
#include "Grid.h"

BeginNameSpace( ONEFLOW )

Unsteady::Unsteady()
{
    sTid = -1;
    data = 0;
    field = 0;
}

Unsteady::~Unsteady()
{
}

void Unsteady::UpdateUnsteady( int sTid )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );

    UsdPara * usdPara = fieldManager->usdPara;

    Grid * grid = Zone::GetGrid();

    MRField * q  = GetFieldPointer< MRField > ( grid, usdPara->flow[ 0 ] );
    MRField * q1 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 1 ] );
    MRField * q2 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 2 ] );

    SetField( q2, q1 );
    SetField( q1, q  );
}


EndNameSpace