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
#include "CpuEulerDomainBackend.h"
#include "EulerDomainStateSync.h"
#include "EulerRungeKuttaCapability.h"
#include "EulerDomainStateRegistry.h"
#include "SimuContext.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "ZoneState.h"
#include "NsCom.h"
#include <stdexcept>

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

namespace
{

struct EulerRungeKuttaStageContext
{
    int stageCount = 0;
};

void RunEulerRungeKuttaStage( int, int stage, void * userData )
{
    auto * context =
        static_cast< EulerRungeKuttaStageContext * >( userData );
    if ( context == nullptr || stage < 0 || stage >= context->stageCount )
    {
        throw std::invalid_argument( "invalid Euler RungeKutta stage callback" );
    }

    ctrl.lhscoef = ctrl.rk_coef[ stage ];
    ONEFLOW::SingleSolverSingleGridTask( "LOAD_RESIDUALS"   );
    ONEFLOW::SingleSolverSingleGridTask( "UPDATE_RESIDUALS" );
    ONEFLOW::SingleSolverSingleGridTask( "CALC_LHS"          );
    ONEFLOW::SingleSolverSingleGridTask( "UPDATE_FLOWFIELD" );
    ONEFLOW::SingleSolverSingleGridTask( "CALC_BOUNDARY"    );
}

EulerRungeKuttaCapabilityRequest CurrentEulerRungeKuttaCapabilityRequest(
    const SimuContext & context )
{
    EulerRungeKuttaCapabilityRequest request;
    request.solverType = SolverState::solverType;
    request.localZoneCount = ZoneState::nLocal;
    request.gridLevel = GridState::gridLevel;
    request.gridCount = GridState::nGrids;
    request.nEquations = nscom.nEqu;
    request.inviscidScheme = nscom.ischeme;
    request.timeIntegral =
        static_cast< EulerRungeKuttaTimeIntegral >( ctrl.time_integral );
    request.hasViscousTerms = nscom.nTModel != 0;
    request.hasSourceTerms = nscom.chemModel != 0;
    request.hasLimiter = ctrl.ilim != 0;
    request.hasInterfaceExchange = ZoneState::nLocal > 1;
    const EulerDomainStateKey key{
        SolverState::solverIndex,
        ZoneState::zid,
        GridState::gridLevel,
        AccelBackendKind::CPU };
    request.backendSupportsAdvance = context.AccelStates().Contains( key );
    return request;
}

}

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
        TimeIntegral::timeIntegral =
            static_cast< TIME_INTEGRAL >( & TimeIntegral::RungeKutta );
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

void TimeIntegral::Relaxation( int nCycles, SimuContext & context )
{
    for ( int iCycle = 0; iCycle < nCycles; ++ iCycle )
    {
        if ( ctrl.time_integral == MULTI_STAGE )
        {
            TimeIntegral::RungeKutta( context );
        }
        else if ( ctrl.time_integral == SIMPLE )
        {
            TimeIntegral::Simple();
        }
        else
        {
            TimeIntegral::Lusgs();
        }
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

void TimeIntegral::RungeKutta( SimuContext & context )
{
    if ( GridState::gridLevel != 0 )
    {
        TimeIntegral::RungeKutta();
        return;
    }

    const EulerRungeKuttaCapabilityDecision decision =
        EvaluateEulerRungeKuttaCapability(
            CurrentEulerRungeKuttaCapabilityRequest( context ) );
    if ( ! decision.enabled || ctrl.rk_coef.size() == 0 )
    {
        // The legacy task sequence remains the numerical fallback until a
        // backend can own conserved-state stages and halo/boundary exchange.
        TimeIntegral::RungeKutta();
        return;
    }

    const EulerDomainStateKey key{
        SolverState::solverIndex,
        ZoneState::zid,
        GridState::gridLevel,
        AccelBackendKind::CPU };
    if ( ! context.AccelStates().Contains( key ) )
    {
        TimeIntegral::RungeKutta();
        return;
    }

    CpuEulerDomainBackend backend;
    EulerDomainState & state = context.AccelStates().Get( key );

    ONEFLOW::SingleSolverSingleGridTask( "LOAD_Q"        );
    ONEFLOW::SingleSolverSingleGridTask( "CALC_TIME_STEP" );

    EulerRungeKuttaStageContext stageContext;
    stageContext.stageCount = static_cast< int >( ctrl.rk_coef.size() );

    EulerDomainRunOptions options;
    options.stageCount = stageContext.stageCount;
    options.stageCallback = & RunEulerRungeKuttaStage;
    options.stageContext = & stageContext;
    backend.Advance( state, 1, options );

    // The existing MRField/task path is authoritative for this first CPU
    // vertical slice; keep the lifecycle cache synchronized at macro-step end.
    UploadCurrentEulerDomainState( context, backend );
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
