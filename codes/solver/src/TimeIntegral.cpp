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

#include "TimeIntegral.h"
#include "Multigrid.h"
#include "CmxTask.h"
#include "CmxTaskNames.h"
#include "Message.h"
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
	else if ( ctrl.time_integral == SIMPLE )
	{
		TimeIntegral::timeIntegral = &TimeIntegral::Simple;
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
    // Resolve source names once per entry; stage/sweep loops use ids only.
    using ONEFLOW::MessageMap;
    const int idLoadQ           = MessageMap::GetMsgId( kLoadQTaskName );
    const int idCalcTimeStep    = MessageMap::GetMsgId( kCalcTimeStepTaskName );
    const int idLoadResiduals   = MessageMap::GetMsgId( kLoadResidualsTaskName );
    const int idUpdateResiduals = MessageMap::GetMsgId( kUpdateResidualsTaskName );
    const int idCalcLhs         = MessageMap::GetMsgId( kCalcLhsTaskName );
    const int idUpdateFlow      = MessageMap::GetMsgId( kUpdateFlowFieldTaskName );
    const int idCalcBoundary    = MessageMap::GetMsgId( kCalcBoundaryTaskName );

    if ( GridState::gridLevel == 0 )
    {
        ONEFLOW::SingleSolverSingleGridTask( idLoadQ );
        ONEFLOW::SingleSolverSingleGridTask( idCalcTimeStep );

        int nStages = ctrl.rk_coef.size();
        for ( int iStage = 0; iStage < nStages; ++ iStage )
        {
            ctrl.lhscoef = ctrl.rk_coef[ iStage ];

            ONEFLOW::SingleSolverSingleGridTask( idLoadResiduals );
            ONEFLOW::SingleSolverSingleGridTask( idUpdateResiduals );
            ONEFLOW::SingleSolverSingleGridTask( idCalcLhs );
            ONEFLOW::SingleSolverSingleGridTask( idUpdateFlow );
            ONEFLOW::SingleSolverSingleGridTask( idCalcBoundary );
        }
    }
    else
    {
        ctrl.lhscoef = 1.0;
        ONEFLOW::SingleSolverSingleGridTask( idLoadQ );
        ONEFLOW::SingleSolverSingleGridTask( idCalcTimeStep );
        ONEFLOW::SingleSolverSingleGridTask( idLoadResiduals );
        ONEFLOW::SingleSolverSingleGridTask( idUpdateResiduals );
        ONEFLOW::SingleSolverSingleGridTask( idCalcLhs );
        ONEFLOW::SingleSolverSingleGridTask( idUpdateFlow );
        ONEFLOW::SingleSolverSingleGridTask( idCalcBoundary );
    }
}

void TimeIntegral::Lusgs()
{
    using ONEFLOW::MessageMap;
    const int idZeroDq          = MessageMap::GetMsgId( kZeroDqFieldTaskName );
    const int idCalcTimeStep    = MessageMap::GetMsgId( kCalcTimeStepTaskName );
    const int idLoadResiduals   = MessageMap::GetMsgId( kLoadResidualsTaskName );
    const int idUpdateResiduals = MessageMap::GetMsgId( kUpdateResidualsTaskName );
    const int idInitLusgs       = MessageMap::GetMsgId( kInitLusgsTaskName );
    const int idLowerSweep      = MessageMap::GetMsgId( kLusgsLowerSweepTaskName );
    const int idExchangeDq      = MessageMap::GetMsgId( kExchangeInterfaceDqTaskName );
    const int idUpperSweep      = MessageMap::GetMsgId( kLusgsUpperSweepTaskName );
    const int idUpdateLusgs     = MessageMap::GetMsgId( kUpdateFlowFieldLusgsTaskName );
    const int idCalcBoundary    = MessageMap::GetMsgId( kCalcBoundaryTaskName );

    ONEFLOW::SingleSolverSingleGridTask( idZeroDq );
    ONEFLOW::SingleSolverSingleGridTask( idCalcTimeStep );
    ONEFLOW::SingleSolverSingleGridTask( idLoadResiduals );
    ONEFLOW::SingleSolverSingleGridTask( idUpdateResiduals );
    ONEFLOW::SingleSolverSingleGridTask( idInitLusgs );

    for ( int iSweep = 0; iSweep < SweepState::nSweeps; ++ iSweep )
    {
        ONEFLOW::SingleSolverSingleGridTask( idLowerSweep );
        ONEFLOW::SingleSolverSingleGridTask( idExchangeDq );
        ONEFLOW::SingleSolverSingleGridTask( idUpperSweep );
    }

    ONEFLOW::SingleSolverSingleGridTask( idUpdateLusgs );
    ONEFLOW::SingleSolverSingleGridTask( idCalcBoundary );
}

void TimeIntegral::Simple()
{
    using ONEFLOW::MessageMap;
    const int idUpdateResiduals = MessageMap::GetMsgId( kUpdateResidualsTaskName );
    const int idSolTurb         = MessageMap::GetMsgId( kSolTurbTaskName );
    const int idSolHeat         = MessageMap::GetMsgId( kSolHeatTaskName );

    ONEFLOW::SingleSolverSingleGridTask( idUpdateResiduals );
    ONEFLOW::SingleSolverSingleGridTask( idSolTurb );
    ONEFLOW::SingleSolverSingleGridTask( idSolHeat );
}

EndNameSpace
