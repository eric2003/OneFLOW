/*---------------------------------------------------------------------------*\\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License either version 3 of the
    License, or (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\\---------------------------------------------------------------------------*/

#include "UUnsteady.h"
#include "UsdData.h"
#include "UsdField.h"
#include "UnsteadyConvergence.h"
#include "TimeIntegration.h"
#include "Iteration.h"
#include "UCom.h"
#include "Com.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

UUnsteady::UUnsteady()
{
}


UUnsteady::~UUnsteady()
{
}

void UUnsteady::UpdateDualTimeStepResidual()
{
    MRField * res =
        field->GetResidual( UsdField::HistoryLevel::Current );

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        ( * res )[ iEqu ][ ug.cId ] =
            dualtimeRes[ iEqu ];
    }
}


void UUnsteady::UpdateDualTimeStepSource()
{
    MRField * res =
        field->GetResidual( UsdField::HistoryLevel::Current );

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        ( * res )[ iEqu ][ ug.cId ] -=
            dualtimeSrc[ iEqu ];
    }
}

void UUnsteady::StoreOldResidual()
{
    //Cxh20140818: first of all, we need to know the residualn1 and residualn2 (the residuals at time n and time n-1);
    //The first step residuals of iteration in two time steps are stored as n-time residuals
    if ( Iteration::innerSteps != 1 ) return;

    MRField * current =
        field->GetResidual( UsdField::HistoryLevel::Current );

    MRField * previous =
        field->GetResidual( UsdField::HistoryLevel::Previous );

    MRField * old =
        field->GetResidual( UsdField::HistoryLevel::Old );

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
        {
            ( * old )[ iEqu ][ cId ] =
                ( * previous )[ iEqu ][ cId ];

            ( * previous )[ iEqu ][ cId ] =
                ( * current )[ iEqu ][ cId ];
        }
    }
}

void UUnsteady::PrepareResidual()
{
    MRField * res =
        field->GetResidual( UsdField::HistoryLevel::Current );

    MRField * res1 =
        field->GetResidual( UsdField::HistoryLevel::Previous );

    MRField * res2 =
        field->GetResidual( UsdField::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        this->res[ iEqu ] =
            ( * res )[ iEqu ][ ug.cId ];

        this->res1[ iEqu ] =
            ( * res1 )[ iEqu ][ ug.cId ];

        this->res2[ iEqu ] =
            ( * res2 )[ iEqu ][ ug.cId ];
    }
}

void UUnsteady::CalcCellDualTimeResidual()
{
    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        dualtimeRes[ iEqu ] = timeIntegration.resc1 * res [ iEqu ] +
                               timeIntegration.resc2 * res1[ iEqu ] +
                               timeIntegration.resc3 * res2[ iEqu ];
    }
}

void UUnsteady::CalcCellDualTimeSrc()
{
    const RealField & q  = this->q;
    const RealField & q1 = this->q1;
    const RealField & q2 = this->q2;

    for ( int iEqu = 0; iEqu < data->nEqu; ++ iEqu )
    {
        Real dualSrc0 = timeIntegration.sc1 * gcom.cvol  * q [ iEqu ];
        Real dualSrc1 = timeIntegration.sc2 * gcom.cvol1 * q1[ iEqu ];
        Real dualSrc2 = timeIntegration.sc3 * gcom.cvol2 * q2[ iEqu ];

        dualtimeSrc[ iEqu ] = dualSrc0 + dualSrc1 + dualSrc2;
    }
}

void UUnsteady::CalcDualTimeResidual()
{
    timeIntegration.CalcResCoef();
    res.resize( data->nEqu );
    res1.resize( data->nEqu );
    res2.resize( data->nEqu );
    dualtimeRes.resize( data->nEqu );

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId = cId;

        this->PrepareResidual();

        this->CalcCellDualTimeResidual();

        this->UpdateDualTimeStepResidual();
    }
}

void UUnsteady::CalcDualTimeSrc()
{
    ug.Init();
    this->StoreOldResidual();

    this->CalcDualTimeResidual();

    timeIntegration.CalcSrcCoeff();
    dualtimeSrc.resize( data->nEqu );

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId = cId;

        ( * this->srcFun )( this );

        this->CalcCellDualTimeSrc();

        this->UpdateDualTimeStepSource();
    }
}

void UUnsteady::CalcUnsteadyCriterion()
{
    data->convergence.Reset();
    res.resize( data->nEqu );

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId = cId;

        ( * this->criFun )( this );

        this->PrepareResidual();

        data->convergence.Accumulate(
            res,
            data->GetQ1(),
            data->GetQ2() );
    }

    data->convergence.Calculate();
}

EndNameSpace
