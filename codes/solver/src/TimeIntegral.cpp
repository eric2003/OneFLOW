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
#include "StageProfiler.h"
#include "SimuContext.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "ZoneState.h"
#include "NsCom.h"
#include "FieldWrap.h"
#include "DataBase.h"
#include "Grid.h"
#include "Zone.h"
#include "Parallel.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

int SweepState::nSweeps = 1;

TIME_INTEGRAL TimeIntegral::timeIntegral;

namespace
{

struct EulerRungeKuttaStageContext
{
    int stageCount = 0;
    SimuContext * context = nullptr;
};

void EmitEulerStageDiagnostic( int stage, const char * taskName )
{
    const char * enabled =
        std::getenv( "ONEFLOW_UNS_STAGE_DIAGNOSTIC" );
    if ( enabled == nullptr || enabled[ 0 ] != '1' ) return;

    Grid * grid = Zone::GetGrid();
    const MRField * q =
        grid == nullptr ? nullptr
                        : GetFieldPointer< MRField >( grid, "q" );
    const MRField * res =
        grid == nullptr ? nullptr
                        : GetFieldPointer< MRField >( grid, "res" );

    auto countNonFinite = []( const MRField * field,
                              int & equation,
                              int & item ) -> std::size_t
    {
        equation = -1;
        item = -1;
        if ( field == nullptr ) return 0;
        std::size_t count = 0;
        for ( int iEqu = 0; iEqu < field->GetNEqu(); ++ iEqu )
        {
            const auto & values = ( * field )[ iEqu ];
            for ( std::size_t i = 0; i < values.size(); ++ i )
            {
                if ( ! std::isfinite( values[ i ] ) )
                {
                    if ( equation < 0 )
                    {
                        equation = iEqu;
                        item = static_cast< int >( i );
                    }
                    ++ count;
                }
            }
        }
        return count;
    };

    int qEquation = -1;
    int qItem = -1;
    int resEquation = -1;
    int resItem = -1;
    const std::size_t qBad = countNonFinite( q, qEquation, qItem );
    const std::size_t resBad =
        countNonFinite( res, resEquation, resItem );

    std::cout << "STAGE_DIAG rank=" << Parallel::GetPid()
              << " zone=" << ZoneState::zid
              << " stage=" << stage
              << " task=" << taskName
              << " q_bad=" << qBad
              << " q_first=" << qEquation << ":" << qItem
              << " res_bad=" << resBad
              << " res_first=" << resEquation << ":" << resItem
              << std::endl;
}

void RunEulerRungeKuttaStage( int, int stage, void * userData )
{
    ScopedStageTimer timer( "rk_update" );
    auto * context =
        static_cast< EulerRungeKuttaStageContext * >( userData );
    if ( context == nullptr || stage < 0 || stage >= context->stageCount )
    {
        throw std::invalid_argument( "invalid Euler RungeKutta stage callback" );
    }

    ctrl.lhscoef = ctrl.rk_coef[ stage ];
    EmitEulerStageDiagnostic( stage, "BEGIN" );
    if ( context->context != nullptr )
    {
        ONEFLOW::SingleSolverSingleGridTask(
            "LOAD_RESIDUALS", *context->context );
        EmitEulerStageDiagnostic( stage, "LOAD_RESIDUALS" );
        ONEFLOW::SingleSolverSingleGridTask(
            "UPDATE_RESIDUALS", *context->context );
        EmitEulerStageDiagnostic( stage, "UPDATE_RESIDUALS" );
        ONEFLOW::SingleSolverSingleGridTask(
            "CALC_LHS", *context->context );
        EmitEulerStageDiagnostic( stage, "CALC_LHS" );
        ONEFLOW::SingleSolverSingleGridTask(
            "UPDATE_FLOWFIELD", *context->context );
        EmitEulerStageDiagnostic( stage, "UPDATE_FLOWFIELD" );
        ONEFLOW::SingleSolverSingleGridTask(
            "CALC_BOUNDARY", *context->context );
        EmitEulerStageDiagnostic( stage, "CALC_BOUNDARY" );
    }
    else
    {
        ONEFLOW::SingleSolverSingleGridTask( "LOAD_RESIDUALS" );
        EmitEulerStageDiagnostic( stage, "LOAD_RESIDUALS" );
        ONEFLOW::SingleSolverSingleGridTask( "UPDATE_RESIDUALS" );
        EmitEulerStageDiagnostic( stage, "UPDATE_RESIDUALS" );
        ONEFLOW::SingleSolverSingleGridTask( "CALC_LHS" );
        EmitEulerStageDiagnostic( stage, "CALC_LHS" );
        ONEFLOW::SingleSolverSingleGridTask( "UPDATE_FLOWFIELD" );
        EmitEulerStageDiagnostic( stage, "UPDATE_FLOWFIELD" );
        ONEFLOW::SingleSolverSingleGridTask( "CALC_BOUNDARY" );
        EmitEulerStageDiagnostic( stage, "CALC_BOUNDARY" );
    }
    // A multi-zone explicit stage must publish updated boundary values before
    // the next stage consumes ghost cells.  The legacy outer-iteration
    // exchange is too late for RK stage semantics.
    if ( ZoneState::nZones > 1 )
    {
        if ( context->context != nullptr )
        {
            ONEFLOW::SingleSolverSingleGridTask(
                "UPDATE_INTERFACE_DATA", *context->context );
        }
        else
        {
            ONEFLOW::SingleSolverSingleGridTask(
                "UPDATE_INTERFACE_DATA" );
        }
        EmitEulerStageDiagnostic( stage, "UPDATE_INTERFACE_DATA" );
    }
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
    request.hasViscousTerms = vis_model.vismodel != 0;
    request.hasSourceTerms = nscom.chemModel != 0;
    request.hasLimiter = ctrl.ilim != 0;
    // A partitioned multi-zone run may assign exactly one zone to each rank.
    // Use the global zone count as well as the local count; localZoneCount
    // alone would incorrectly enable the resident path for cross-rank halos.
    request.hasInterfaceExchange =
        ZoneState::nZones > 1 || ZoneState::nLocal > 1;
    const EulerDomainStateKey key{
        SolverState::solverIndex,
        ZoneState::zid,
        GridState::gridLevel,
        AccelBackendKind::CPU };
    request.backendSupportsAdvance = context.AccelStates().Contains( key );
    return request;
}

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
            EmitEulerStageDiagnostic( iStage, "LOAD_RESIDUALS" );
            ONEFLOW::SingleSolverSingleGridTask( idUpdateResiduals );
            EmitEulerStageDiagnostic( iStage, "UPDATE_RESIDUALS" );
            ONEFLOW::SingleSolverSingleGridTask( idCalcLhs );
            EmitEulerStageDiagnostic( iStage, "CALC_LHS" );
            ONEFLOW::SingleSolverSingleGridTask( idUpdateFlow );
            EmitEulerStageDiagnostic( iStage, "UPDATE_FLOWFIELD" );
            ONEFLOW::SingleSolverSingleGridTask( idCalcBoundary );
            EmitEulerStageDiagnostic( iStage, "CALC_BOUNDARY" );
            if ( ZoneState::nZones > 1 )
            {
                ONEFLOW::SingleSolverSingleGridTask(
                    kUpdateInterfaceDataTaskName );
                EmitEulerStageDiagnostic(
                    iStage, "UPDATE_INTERFACE_DATA" );
            }
        }
    }
    else
    {
        ctrl.lhscoef = 1.0;
        ONEFLOW::SingleSolverSingleGridTask( idLoadQ );
        ONEFLOW::SingleSolverSingleGridTask( idCalcTimeStep );
        ONEFLOW::SingleSolverSingleGridTask( idLoadResiduals );
        EmitEulerStageDiagnostic( 0, "LOAD_RESIDUALS" );
        ONEFLOW::SingleSolverSingleGridTask( idUpdateResiduals );
        EmitEulerStageDiagnostic( 0, "UPDATE_RESIDUALS" );
        ONEFLOW::SingleSolverSingleGridTask( idCalcLhs );
        EmitEulerStageDiagnostic( 0, "CALC_LHS" );
        ONEFLOW::SingleSolverSingleGridTask( idUpdateFlow );
        EmitEulerStageDiagnostic( 0, "UPDATE_FLOWFIELD" );
        ONEFLOW::SingleSolverSingleGridTask( idCalcBoundary );
        EmitEulerStageDiagnostic( 0, "CALC_BOUNDARY" );
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

    ONEFLOW::SingleSolverSingleGridTask( "LOAD_Q", context );
    ONEFLOW::SingleSolverSingleGridTask( "CALC_TIME_STEP", context );

    EulerRungeKuttaStageContext stageContext;
    stageContext.stageCount = static_cast< int >( ctrl.rk_coef.size() );
    stageContext.context = &context;

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
