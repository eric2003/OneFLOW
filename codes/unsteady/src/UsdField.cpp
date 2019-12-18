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

#include "UsdField.h"
#include "UsdPara.h"
#include "FieldImp.h"
#include "FieldWrap.h"
#include "DataBase.h"
#include "Zone.h"
#include "UnsGrid.h"

BeginNameSpace( ONEFLOW )


UsdField::UsdField()
{
}

UsdField::~UsdField()
{
}

void UsdField::Init()
{
    ;
}

void UsdField::InitBasic( int sTid )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );
    UsdPara * usdPara = fieldManager->usdPara;
    q  = GetFieldPointer< MRField > ( grid, usdPara->flow[ 0 ] );
    q1 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 1 ] );
    q2 = GetFieldPointer< MRField > ( grid, usdPara->flow[ 2 ] );

    res  = GetFieldPointer< MRField > ( grid, usdPara->residual[ 0 ] );
    res1 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 1 ] );
    res2 = GetFieldPointer< MRField > ( grid, usdPara->residual[ 2 ] );
}


EndNameSpace