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

#include "NsSolverImp.h"
#include "SolverImp.h"
#include "UTimeStep.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "GridState.h"
#include "UNsVisFlux.h"
#include "NsRhs.h"
#include "UCom.h"
#include "CmxTask.h"
#include "CmxTaskNames.h"
#include "Iteration.h"
#include "SolverDef.h"
#include "Ctrl.h"
#include "TimeIntegral.h"
#include "LaminarPlate.h"
#include "TurbPlate.h"
#include "SolverRegister.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "StageProfiler.h"
#include "AccelRuntime.h"
#include "SolverState.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "UNsCom.h"
// The policy header is header-only and CPU-safe (it returns false when the
// HIP backend is not built/selected). The boundary gate calls it in both
// builds, so it must not be guarded by ONEFLOW_ENABLE_HIP.
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <cstdlib>
#include <iostream>
#include <string>
#include <iomanip>
#include <memory>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

namespace
{

bool HipBoundaryExplicitlyRequested()
{
    return HipProductionFeatureExplicitlyEnabled(
        "ONEFLOW_ENABLE_UNS_HIP_BOUNDARY" );
}

bool HipBackendSelected()
{
#ifdef ONEFLOW_ENABLE_HIP
    AccelRuntime & runtime = AccelRuntime::Instance();
    return runtime.IsInitialized()
        && runtime.Backend().Kind() == AccelBackendKind::HIP;
#else
    return false;
#endif
}

bool HipBoundarySupported()
{
#ifndef ONEFLOW_ENABLE_HIP
    return false;
#else
    // One-shot rejection reason, enabled by the same diagnostic switch as the
    // gradient/reconstruction decision line. A silent false would otherwise be
    // indistinguishable from "the user did not ask for the device path".
    auto reject = []( const char * reason )
    {
        if ( std::getenv( "ONEFLOW_UNS_HIP_STAGE_DIAGNOSTIC" ) == nullptr )
            return false;
        static std::string reported;
        if ( reported.find( reason ) == std::string::npos )
        {
            reported += reason;
            reported += ';';
            std::cout << "HIP_BOUNDARY_DIAGNOSTIC unsupported=" << reason
                      << std::endl;
        }
        return false;
    };

    AccelRuntime & runtime = AccelRuntime::Instance();
    if ( ! runtime.IsInitialized()
         || runtime.Backend().Kind() != AccelBackendKind::HIP )
        return reject( "runtime" );
    // Multi-zone is allowed: the device boundary kernel skips the BC types
    // that an external halo exchange owns (INTERFACE -1, PERIODIC 99,
    // OVERSET 1000), so interface ghosts stay owned by the interface contract
    // while physical boundary ghosts are evaluated on the device.
    if ( SolverState::solverType != NS_SOLVER )
        return reject( "solver" );
    if ( nscom.nEqu != 5 ) return reject( "equations" );
    if ( vis_model.vismodel != 0 ) return reject( "viscous" );
    if ( ctrl.ilim != 0 ) return reject( "limiter" );
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr ) return reject( "grid" );
    // NsCalcBoundary runs per local zone; in a partitioned run the cached
    // unstructured view may still point at the previous zone, so re-bind it
    // before validating the boundary set.
    ug.Init();
    if ( ug.nBFaces <= 0 ) return reject( "no_boundary_faces" );
    if ( ug.nTCell - ug.nCells != ug.nBFaces )
        return reject( "ghost_extent" );
    if ( ug.bcRecord == nullptr ) return reject( "bc_record" );
    for ( int face = 0; face < ug.nBFaces; ++ face )
    {
        const int type = ug.bcRecord->bcType[ face ];
        if ( type < 0 || type == BC::PERIODIC || type == BC::OVERSET
             || type == BC::EXTRAPOLATION || type == BC::SOLID_SURFACE
             || type == BC::SYMMETRY || type == BC::FARFIELD
             || type == BC::INFLOW || type == BC::OUTFLOW
             || type == BC::POLE || type == BC::POLE1
             || type == BC::POLE2 || type == BC::POLE3 )
            continue;
        static std::string bcReason;
        if ( bcReason.empty() )
            bcReason = "bc_type_" + std::to_string( type );
        return reject( bcReason.c_str() );
    }
    return true;
#endif
}

#ifdef ONEFLOW_ENABLE_HIP
void SyncHipBoundaryGhostFromHost()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr ) return;
    MRField * q = GetFieldPointer< MRField >( grid, "q" );
    if ( q == nullptr ) return;
    const int nEquations = static_cast< int >( q->GetNEqu() );
    if ( nEquations < 3 || nEquations > 5 ) return;
    const int nGhostCells = grid->nBFaces;
    if ( nGhostCells <= 0 ) return;
    const Real * components[ 5 ] = {};
    for ( int equation = 0; equation < nEquations; ++ equation )
        components[ equation ] = ( * q )[ equation ].data();

    HipFluxBackend & backend = HipFluxBackend::Shared();
    // The host boundary model owns the physical ghost update.  Publish the
    // mutation and, when the resident ghost epoch is stale, refresh only the
    // ghost region of the state-owned device q.  This keeps the device ghost
    // bit-identical to the host ghost without re-uploading the whole field.
    backend.NotifyHostGhostStateModified( grid );
    backend.SyncResidentGhostStateFromHost(
        grid, grid->nCells, nGhostCells, nEquations, components );
}

// True when this step runs the device-authoritative chain: the device owns the
// newest q and the host mirror is only refreshed on demand.
bool HipDeviceAuthoritativeHostMirror()
{
    return HipDeviceAuthoritativeStateEnabled()
        && HipProductionBackendSelected()
        && ctrl.startStrategy != 2;
}
#endif

void CalcHipBoundary()
{
#ifdef ONEFLOW_ENABLE_HIP
    UnsGrid * grid = Zone::GetUnsGrid();
    // Only skip the per-stage ghost mirror once the device really owns a newer
    // field than the host.  During initialization (and before the first
    // device advance) the host mirror is the only complete source, and the
    // first full upload of the bound state must carry valid ghost cells.
    const bool deviceAuthoritative =
        HipDeviceAuthoritativeHostMirror()
        && HipFluxBackend::Shared().HostMirrorStale( grid )
        && ! HipGhostQMirrorForced();
    BoundaryStateView view;
    {
    ScopedStageTimer prepTimer( "bnd_prep" );
    ug.Init();
    unsf.Init();
    view.nCells = ug.nCells;
    view.nGhostCells = ug.nTCell - ug.nCells;
    view.nBoundaryFaces = ug.nBFaces;
    view.nEquations = nscom.nEqu;
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    view.boundaryType = ug.bcRecord->bcType.data();
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.meshVelocityNormal = ug.vfn->data();
    view.inflow = nscom.inflow.data();
    view.gamma = nscom.gama_ref;
    view.faceOuterNormal = nscom.faceOuterNormal;
    view.cacheKey = grid;
    view.preferResidentState = true;
    // The device owns the newest q: the ghost cells stay on the device while
    // the wall state keeps its per-stage host mirror (the reconstruction
    // consumes that host copy).  The transfer can be forced back on for
    // bisection/diagnostics.
    view.deviceAuthoritativeHostMirror = deviceAuthoritative;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
    {
        view.q[ equation ] = ( * unsf.q )[ equation ].data();
        if ( unsf.bc_q != nullptr )
            view.bcQ[ equation ] = ( * unsf.bc_q )[ equation ].data();
    }
    }
    {
    ScopedStageTimer backendTimer( "bnd_backend" );
    HipFluxBackend::Shared().UpdateBoundaryState( view );
    }
#else
    throw std::runtime_error(
        "HIP boundary update was requested without HIP support" );
#endif
}

}

void RegisterNsFunc()
{
    REGISTER_DATA_CLASS( NsInitFinal );
    REGISTER_DATA_CLASS( NsVisual );
    REGISTER_DATA_CLASS( NsCalcTimeStep );
    REGISTER_DATA_CLASS( NsUpdateResiduals );
    REGISTER_DATA_CLASS( NsImplicitMethod );
    REGISTER_DATA_CLASS( NsPostprocess );
    REGISTER_DATA_CLASS( NsFinalPostprocess );
    REGISTER_DATA_CLASS( NsInitSolver );
    REGISTER_DATA_CLASS( NsCalcBoundary );
    REGISTER_DATA_CLASS( DumpHeatFluxCoeff );
    SetPlateTask();
    SetTurbPlateTask();
}

void NsInitFinal( StringField & data )
{
    //NsCalcGamaT( F_INNER );
    //CalcLaminarViscosity( F_INNER );
    //NsCalcBc();
    //NsCalcGamaT( F_GHOST );
    //CalcLaminarViscosity( F_GHOST );

    NsCalcBoundary( data );

    Grid & grid = Zone::GetGridReference();

    if ( Zone::GetCGrid( &grid ) )
    {
        //RestrictAllQ( NS_SOLVER, FLOW_FIELD_INDEX );

        GridState::gridLevel += 1;

        NsCalcBoundary( data );

        GridState::gridLevel -= 1;
    }
}

void NsVisual( StringField & data )
{
    ;
}

void NsCalcBoundary( StringField & data )
{
    if ( vis_model.vismodel > 0 && nscom.nTModel <= 0 )
    {
        throw std::runtime_error(
            "viscous or turbulent NS boundary conditions require "
            "nTModel > 0" );
    }
    const bool explicitRequest = HipBoundaryExplicitlyRequested();
    // The boundary stage follows the backend policy; the capability check
    // still decides, and an explicit request that is unsupported fails fast
    // instead of silently running the host model.
    const bool requested =
        HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_BOUNDARY" )
        && HipBackendSelected() && HipBoundarySupported();
    if ( explicitRequest && ! requested )
    {
        throw std::runtime_error(
            "ONEFLOW_ENABLE_UNS_HIP_BOUNDARY=1 is unsupported for the "
            "current solver, boundary set, or runtime" );
    }
#ifdef ONEFLOW_ENABLE_HIP
    // With a device-authoritative q the host mirror is stale, so the host
    // gamma/temperature fields would be rebuilt from old values once per
    // stage.  The device path uses the solver gamma constant, and the derived
    // host fields are rebuilt when a host reader refreshes the mirror.
    const bool skipHostDerivedFields =
        requested && HipDeviceAuthoritativeHostMirror()
        && HipFluxBackend::Shared().HostMirrorStale( Zone::GetUnsGrid() );
#else
    const bool skipHostDerivedFields = false;
#endif
    {
        ScopedStageTimer timer( "boundary_gamma_inner" );
        if ( ! skipHostDerivedFields )
        {
#ifdef ONEFLOW_ENABLE_HIP
            // The host model (and the derived gamma/temperature fields) reads
            // the host q mirror; make sure it is current first.
            if ( HipBackendSelected() ) NsRefreshHostStateMirror();
#endif
            NsCalcGamaT( F_INNER );
        }
    }
    if ( vis_model.vismodel != 0 )
    {
        ScopedStageTimer timer( "boundary_viscosity_inner" );
        CalcLaminarViscosity( F_INNER );
    }
    if ( requested )
    {
        ScopedStageTimer timer( "boundary_device" );
        CalcHipBoundary();
    }
    else
    {
        ScopedStageTimer timer( "boundary_bc" );
        NsCalcBc();
#ifdef ONEFLOW_ENABLE_HIP
        if ( HipBackendSelected() ) SyncHipBoundaryGhostFromHost();
#endif
    }
    {
        ScopedStageTimer timer( "boundary_gamma_ghost" );
        if ( ! skipHostDerivedFields ) NsCalcGamaT( F_GHOST );
    }
    if ( vis_model.vismodel != 0 )
    {
        ScopedStageTimer timer( "boundary_viscosity_ghost" );
        CalcLaminarViscosity( F_GHOST );
    }
}

void NsCalcTimeStep( StringField & data )
{
    auto uTimeStep = std::make_unique<UTimeStep>();
    uTimeStep->CalcTimeStep();
}

void NsUpdateResiduals( StringField & data )
{
    auto rhs = std::make_unique<NsRhs>();
    rhs->UpdateResiduals();
}

void NsImplicitMethod( StringField & data )
{
    ;
}

namespace {

// The explicit RK path exchanges the interface at the end of every stage, so
// q and its resident device ghosts are already current when the outer
// iteration reaches the post-process hook.  Repeating the upload/exchange/
// download round trip here costs a full interface pass per step.  Opt in with
// ONEFLOW_SKIP_POST_INTERFACE=1 until the trace and residual gates have run.
bool SkipRedundantPostInterface()
{
    const char * value = std::getenv( "ONEFLOW_SKIP_POST_INTERFACE" );
    if ( value == nullptr || value[ 0 ] != '1' ) return false;
    return ctrl.time_integral == MULTI_STAGE
        && ZoneState::nZones > 1;
}

}

void NsPostprocess( StringField & data )
{
    //After every Iter, the first thing to consider is communication.
    if ( ! SkipRedundantPostInterface() )
    {
        CommInterfaceData();
    }

    //The solution and output of residuals need to be judged logically.
    if ( Iteration::ResOk() )
    {
        ONEFLOW::AddCmdToList( kDumpResidualTaskName );
    }

    //The solution and output of aerodynamic force need to be judged logically.
    if ( Iteration::ForceOk() )
    {
        ONEFLOW::AddCmdToList( kDumpAerodynamicTaskName );
    }

    if ( ! Iteration::InnerOk() ) return;

    ONEFLOW::AddCmdToList( kUpdateUnsteadyFlowTaskName );

    if ( Iteration::outerSteps % Iteration::nVisualSave == 0 )
    {
        ONEFLOW::AddCmdToList( kVisualizationTaskName      );
        ONEFLOW::AddCmdToList( kDumpPressureCoeffTaskName  );
        ONEFLOW::AddCmdToList( kDumpHeatfluxCoeffTaskName  );
        if ( ctrl.idump == 1 )
        {
            ONEFLOW::AddCmdToList( kDumpLaminarPlateTaskName );
        }
        else if ( ctrl.idump == 2 )
        {
            ONEFLOW::AddCmdToList( kDumpTurbPlateTaskName );
        }
    }

    if ( Iteration::outerSteps % Iteration::nFieldSave == 0 )
    {
        ONEFLOW::AddCmdToList( kDumpRestartTaskName );
    }

}

void DumpHeatFluxCoeff( StringField & data )
{
    ;
}

void NsFinalPostprocess( StringField & data )
{
    ONEFLOW::AddCmdToList( kDumpRestartTaskName       );
    ONEFLOW::AddCmdToList( kDumpAerodynamicTaskName   );
    ONEFLOW::AddCmdToList( kDumpPressureCoeffTaskName );
    ONEFLOW::AddCmdToList( kDumpHeatfluxCoeffTaskName );
    ONEFLOW::AddCmdToList( kVisualizationTaskName     );
}

void NsInitSolver( StringField & data )
{
    ONEFLOW::CommInterfaceData();
}

SolverRegData nsReg;
SolverRegData * GetNsReg()
{
    nsReg.solverType = NS_SOLVER;
    nsReg.func = & RegisterNsFunc;
    nsReg.solverName = "ns";
    nsReg.baseKind = SOLVER_BASED;
    nsReg.dataFlag = WITH_DATA;
    return & nsReg;
}

REGISTER_REG_DATA( GetNsReg );

EndNameSpace
