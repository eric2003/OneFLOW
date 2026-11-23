/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "NsResidentMacroStep.h"

#include "NsSolverImp.h"
#include "ULhs.h"
#include "UNsUpdate.h"
#include "UTimeStep.h"
#include "Ctrl.h"
#include "GridState.h"
#include "Multigrid.h"
#include "HipProductionPolicy.h"
#include "HipFluxBackend.h"
#include "ResidualBasePolicy.h"
#include "StageProfiler.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "ZoneState.h"
#include "CmxTask.h"
#include "CmxTaskNames.h"

#include <cstdlib>

BeginNameSpace( ONEFLOW )

namespace
{

bool ResidentMacroStepRequested()
{
    // Default on for the configuration this driver covers.  The task path stays
    // one environment variable away (ONEFLOW_HIP_RESIDENT_MACROSTEP=0) for
    // bisection, and remains the path for every configuration outside the
    // device-owned explicit chain.
    const char * value =
        std::getenv( "ONEFLOW_HIP_RESIDENT_MACROSTEP" );
    if ( value != nullptr ) return value[ 0 ] != '0';
    return true;
}

}

bool NsResidentMacroStepSupported()
{
    if ( ! ResidentMacroStepRequested() ) return false;
    // One macro-step owns the device state for the whole RK loop.  Coarse grids
    // and multigrid plumbing stay on the task path.
    if ( GridState::gridLevel != 0 || MG::nMulti > 1 ) return false;
    if ( SolverState::solverType != NS_SOLVER ) return false;
    return DeviceOwnsExplicitRkBookkeeping();
}

void NsResidentMacroStep()
{
    ScopedHostTimer hostTimer( "resident_macrostep" );
    StringField unused;

    // CALC_TIME_STEP: local spectral radius and time step for this step.
    // The macro step bypasses the task path, so this work had no profiler
    // category at all: every breakdown file reported calc_time_step = 0 while
    // the per-step base was several milliseconds.  Time it here so the A1.7
    // base accounting can name the segment instead of inferring it.
    {
        ScopedStageTimer timer( "macro_time_step" );
        ScopedHostTimer hostTimer( "macro_time_step" );
        UTimeStep timeStep;
        timeStep.CalcTimeStep();
    }

    const int nStages = static_cast< int >( ctrl.rk_coef.size() );
    for ( int stage = 0; stage < nStages; ++ stage )
    {
        ctrl.rkStage = stage;
        ctrl.lhscoef = ctrl.rk_coef[ stage ];

        // UPDATE_RESIDUALS: gradient, reconstruction, fused face flux, RHS.
        {
            ScopedStageTimer timer( "macro_update_residuals" );
            ScopedHostTimer hostTimer( "macro_update_residuals" );
            NsUpdateResiduals( unused );
        }
        // CALC_LHS: scale the resident residual by dt/vol*rkcoef.
        {
            ScopedStageTimer timer( "macro_calc_lhs" );
            ScopedHostTimer hostTimer( "macro_calc_lhs" );
            ULhs lhs;
            lhs.CalcLHS( NS_SOLVER );
        }
        // UPDATE_FLOWFIELD: primitive state update on the device.
        {
            ScopedStageTimer timer( "macro_update_flowfield" );
            ScopedHostTimer hostTimer( "macro_update_flowfield" );
            UNsUpdate update;
            update.UpdateFlowField( NS_SOLVER );
        }
        // CALC_BOUNDARY: local boundary models on the resident state.
        {
            ScopedStageTimer timer( "macro_calc_boundary" );
            ScopedHostTimer hostTimer( "macro_calc_boundary" );
            NsCalcBoundary( unused );
        }
        // UPDATE_INTERFACE_DATA: the explicit stage must publish updated
        // boundary values before the next stage reads its ghosts.  Still one
        // task dispatch per stage (not four), because the exchange itself is
        // owned by the interface task and its lag must stay exactly as in the
        // validated task path.
        if ( ZoneState::nZones > 1 && ! DiagNoCommunication() )
        {
            ScopedStageTimer timer( "macro_interface" );
            ScopedHostTimer hostTimer( "macro_interface" );
            SingleSolverSingleGridTask( kUpdateInterfaceDataTaskName );
        }
    }
    ctrl.rkStage = -1;

#ifdef ONEFLOW_ENABLE_HIP
    // With ONEFLOW_HIP_DEFER_STATUS=1 the non-physical/non-finite flags stay on
    // the device while the stages run.  The check is requested here and
    // consumed by the first host-staged interface pack, which drains the queue
    // for its own copy; a configuration without interfaces falls back to the
    // blocking read so the per-step cadence is preserved either way.
    if ( HipStatusCheckDeferred() )
    {
        if ( ZoneState::nZones > 1 && ! DiagNoCommunication() )
        {
            HipFluxBackend::Shared().RequestResidentStatusCheck();
        }
        else
        {
            HipFluxBackend::Shared().CheckResidentStatus();
        }
    }
#endif
}

EndNameSpace
