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
#include "TimeIntegration.h"
#include "Iteration.h"
#include "UCom.h"
#include "Com.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

UUnsteady::UUnsteady()
{
    timeIntegration.Init();
}


int UUnsteady::GetEquationCount() const
{
    return nEqu;
}

RealField & UUnsteady::GetPrimitive(
    UnsteadyFieldView::HistoryLevel level )
{
    switch ( level )
    {
    case UnsteadyFieldView::HistoryLevel::Current:
        return prim;

    case UnsteadyFieldView::HistoryLevel::Previous:
        return prim1;

    case UnsteadyFieldView::HistoryLevel::Old:
        return prim2;
    }

    return prim;
}

RealField & UUnsteady::GetConservative(
    UnsteadyFieldView::HistoryLevel level )
{
    switch ( level )
    {
    case UnsteadyFieldView::HistoryLevel::Current:
        return q;

    case UnsteadyFieldView::HistoryLevel::Previous:
        return q1;

    case UnsteadyFieldView::HistoryLevel::Old:
        return q2;
    }

    return q;
}

void UUnsteady::SetSourceFunction( USDFunc function )
{
    srcFun = function;
}

void UUnsteady::SetCriterionFunction( USDFunc function )
{
    criFun = function;
}

void UUnsteady::SetEquationCount( int equationCount )
{
    nEqu = equationCount;
    prim.resize( nEqu );
    prim1.resize( nEqu );
    prim2.resize( nEqu );
}

void UUnsteady::UpdateDualTimeStepResidual()
{
    MRField * res =
        GetResidual( UnsteadyFieldView::HistoryLevel::Current );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        ( * res )[ iEqu ][ ug.cId ] =
            dualtimeRes[ iEqu ];
    }
}


void UUnsteady::UpdateDualTimeStepSource()
{
    MRField * res =
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Current );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
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
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Current );

    MRField * previous =
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Previous );

    MRField * old =
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Old );

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
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
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Current );

    MRField * res1 =
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Previous );

    MRField * res2 =
        field.GetResidual( UnsteadyFieldView::HistoryLevel::Old );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
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
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dualtimeRes[ iEqu ] = timeIntegration.resc1 * res [ iEqu ] +
                               timeIntegration.resc2 * res1[ iEqu ] +
                               timeIntegration.resc3 * res2[ iEqu ];
    }
}

void UUnsteady::CalcCellDualTimeSrc()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
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
    res.resize( nEqu );
    res1.resize( nEqu );
    res2.resize( nEqu );
    dualtimeRes.resize( nEqu );

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
    dualtimeSrc.resize( nEqu );
    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

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
    convergence.Init( nEqu );
    convergence.Reset();
    res.resize( nEqu );
    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId = cId;

        ( * this->criFun )( this );

        this->PrepareResidual();

        convergence.Accumulate(
            res,
            q1,
            q2 );
    }

    convergence.Calculate();
}

EndNameSpace
