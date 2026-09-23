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

void UsdField::InitBasic( int solverType )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    FieldManager * fieldManager =
        FieldFactory::GetFieldManager( solverType );

    UsdPara * usdPara =
        &fieldManager->GetUsdPara();

    this->flow.resize( usdPara->flow.size() );

    for ( std::size_t i = 0; i < usdPara->flow.size(); ++ i )
    {
        this->flow[ i ] =
            GetFieldPointer< MRField >(
                grid,
                usdPara->flow[ i ] );
    }

    this->residual.resize( usdPara->residual.size() );

    for ( std::size_t i = 0; i < usdPara->residual.size(); ++ i )
    {
        this->residual[ i ] =
            GetFieldPointer< MRField >(
                grid,
                usdPara->residual[ i ] );
    }

    q  = this->flow[ 0 ];
    q1 = this->flow[ 1 ];
    q2 = this->flow[ 2 ];

    res  = this->residual[ 0 ];
    res1 = this->residual[ 1 ];
    res2 = this->residual[ 2 ];
}


EndNameSpace
