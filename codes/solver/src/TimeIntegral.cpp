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

#include "TimeIntegral.h"
#include "Multigrid.h"
#include "CmxTask.h"
#include "GridState.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

int SweepState::nSweeps = 1;

SweepState::SweepState()
{
    ;
}

SweepState::~SweepState()
{
    ;
}

TIME_INTEGRAL TimeIntegral::timeIntegral;

TimeIntegral::TimeIntegral()
{
    ;
}

TimeIntegral::~TimeIntegral()
{
    ;
}

void TimeIntegral::Init()
{
    if ( ctrl.time_integral == MULTI_STAGE )
    {
        TimeIntegral::timeIntegral = & TimeIntegral::RungeKutta;
    }
    else
    {
        TimeIntegral::timeIntegral = & TimeIntegral::Lusgs;
    }
    
}

void TimeIntegral::Relaxation( int nCycles )
{
    TimeIntegral::Init();
    for ( int iCycle = 0; iCycle < nCycles; ++ iCycle )
    {
        TimeIntegral::timeIntegral();
    }
}

void TimeIntegral::RungeKutta()
{
    if ( GridState::gridLevel == 0 )
    {
        ONEFLOW::SsSgTask( "LOAD_Q"        );
        ONEFLOW::SsSgTask( "CALC_TIME_STEP" );

        int nStages = ctrl.rk_coef.size();
        for ( int iStage = 0; iStage < nStages; ++ iStage )
        {
            ctrl.lhscoef = ctrl.rk_coef[ iStage ];

            ONEFLOW::SsSgTask( "LOAD_RESIDUALS"   );
            ONEFLOW::SsSgTask( "UPDATE_RESIDUALS" );
            ONEFLOW::SsSgTask( "CALC_LHS"         );
            ONEFLOW::SsSgTask( "UPDATE_FLOWFIELD" );
            ONEFLOW::SsSgTask( "CALC_BOUNDARY"    );
        }
    }
    else
    {
        ctrl.lhscoef = 1.0;
        ONEFLOW::SsSgTask( "LOAD_Q"           );
        ONEFLOW::SsSgTask( "CALC_TIME_STEP"    );
        ONEFLOW::SsSgTask( "LOAD_RESIDUALS"   );
        ONEFLOW::SsSgTask( "UPDATE_RESIDUALS" );
        ONEFLOW::SsSgTask( "CALC_LHS"         );
        ONEFLOW::SsSgTask( "UPDATE_FLOWFIELD" );
        ONEFLOW::SsSgTask( "CALC_BOUNDARY"    );
    }
}

void TimeIntegral::Lusgs()
{
    ONEFLOW::SsSgTask( "ZERO_DQ_FIELD"    );
    ONEFLOW::SsSgTask( "CALC_TIME_STEP"   );
    ONEFLOW::SsSgTask( "LOAD_RESIDUALS"   );
    ONEFLOW::SsSgTask( "UPDATE_RESIDUALS" );
    ONEFLOW::SsSgTask( "INIT_LUSGS"       );

    for ( int iSweep = 0; iSweep < SweepState::nSweeps; ++ iSweep )
    {
        ONEFLOW::SsSgTask( "LUSGS_LOWER_SWEEP"     );
        ONEFLOW::SsSgTask( "EXCHANGE_INTERFACE_DQ" );
        ONEFLOW::SsSgTask( "LUSGS_UPPER_SWEEP"     );
    }

    ONEFLOW::SsSgTask( "UPDATE_FLOWFIELD_LUSGS" );
    ONEFLOW::SsSgTask( "CALC_BOUNDARY"          );
}

EndNameSpace