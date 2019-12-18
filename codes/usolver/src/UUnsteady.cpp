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

#include "UUnsteady.h"
#include "UsdData.h"
#include "UsdField.h"
#include "Iteration.h"
#include "UCom.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

UUnsteady::UUnsteady()
{
}

UUnsteady::~UUnsteady()
{
}

void UUnsteady::UpdateDualTimeStepResidual()
{
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        ( * field->res )[ iEqu ][ ug.cId ] = data->dualtimeRes[ iEqu ];
    }
}

void UUnsteady::UpdateDualTimeStepSource()
{
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        ( * field->res )[ iEqu ][ ug.cId ] -= data->dualtimeSrc[ iEqu ];
    }
}

void UUnsteady::StoreOldResidual()
{
    //cxh20140818:����Ҫ֪��ResidualN1��ResidualN2��nʱ�̺�n-1ʱ�̵Ĳв��
    //˫ʱ�䲽�ڵ�����һ���в�洢Ϊnʱ�̲в�
    if ( Iteration::innerSteps != 1 ) return;

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
        {
            ( * field->res2 )[ iEqu ][ cId ] = ( * field->res1 )[ iEqu ][ cId ];
            ( * field->res1 )[ iEqu ][ cId ] = ( * field->res  )[ iEqu ][ cId ];
        }
    }
}

void UUnsteady::PrepareResidual()
{
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        data->res [ iEqu ] = ( * field->res  )[ iEqu ][ ug.cId ];
        data->res1[ iEqu ] = ( * field->res1 )[ iEqu ][ ug.cId ];
        data->res2[ iEqu ] = ( * field->res2 )[ iEqu ][ ug.cId ];
    }
}

void UUnsteady::CmpDualTimeResidual()
{
    data->CmpResCoef();

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        this->PrepareResidual();

        data->CmpCellDualTimeResidual();

        this->UpdateDualTimeStepResidual();
    }
}

void UUnsteady::CmpDualTimeSrc()
{
    ug.Init();
    this->StoreOldResidual();

    this->CmpDualTimeResidual();

    data->CmpSrcCoeff();

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        ( * this->srcFun )( this );

        data->CmpCellDualTimeSrc();

        this->UpdateDualTimeStepSource();
    }
}

void UUnsteady::CmpUnsteadyCriterion()
{
    data->ZeroData();

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        ( * this->criFun )( this );
        
        data->CmpCellUnsteadyCri();
    }

    data->CmpCvg();
}

EndNameSpace