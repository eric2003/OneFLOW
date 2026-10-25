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

#include "UTimeStep.h"
#include "NsCom.h"
#include "UNsCom.h"
#include "UCom.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "NsCtrl.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "HXMath.h"
#include "Fatal.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "ZoneSharedFields.h"
#include "SolverDef.h"
#include "FieldBase.h"
#include "Iteration.h"
#include "Ctrl.h"
#include "GridState.h"
#include "ZoneState.h"
#include "SolverState.h"
#include "TimeIntegral.h"
#include "AccelRuntime.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <cstdlib>
#include <iostream>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

namespace
{

bool HipTimeStepRequested()
{
    return HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_TIME_STEP" );
}

bool EnvironmentEnabled( const char * name )
{
    const char * value = std::getenv( name );
    return value != nullptr && value[ 0 ] == '1';
}

bool HipTimeStepSupported()
{
#ifndef ONEFLOW_ENABLE_HIP
    return false;
#else
    AccelRuntime & runtime = AccelRuntime::Instance();
    return runtime.IsInitialized()
        && runtime.Backend().Kind() == AccelBackendKind::HIP
        && SolverState::solverType == NS_SOLVER
        // Multi-zone is supported: the spectral radius reads the resident q
        // ghost, which the resident-state contract keeps equal to the host
        // mirror, and the timestep stays device-resident for the following
        // residual scaling.
        && GridState::gridLevel == 0
        && GridState::nGrids == 1
        && nscom.nEqu == 5
        && nscom.chemModel == 0
        && nscom.nTModel == 1
        && vis_model.vismodel == 0
        && nscom.timestepModel == 0
        && ctrl.time_integral == MULTI_STAGE
        && ctrl.idualtime == 0
        && ug.nCells > 0
        && ug.nFaces > 0
        && ug.nTCell - ug.nCells == ug.nBFaces
        && HipProductionFeatureEnabled(
            "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE" )
        && HipProductionFeatureEnabled(
            "ONEFLOW_ENABLE_UNS_HIP_BOUNDARY" );
#endif
}

void CalcHipLocalInviscidTimeStep()
{
#ifdef ONEFLOW_ENABLE_HIP
    UnsGrid * grid = Zone::GetUnsGrid();
    CellTimeStepView view;
    view.nCells = ug.nCells;
    view.nGhostCells = ug.nTCell - ug.nCells;
    view.nFaces = ug.nFaces;
    view.nBoundaryFaces = ug.nBFaces;
    view.nEquations = nscom.nEqu;
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.meshVelocityNormal = ug.vfn->data();
    view.faceArea = ug.farea->data();
    view.cellVolume = ug.cvol->data();
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    view.inviscidSpectrum = ( * unsf.invsr )[ 0 ].data();
    view.timeStep = ( * unsf.timestep )[ 0 ].data();
    view.gamma = nscom.gama_ref;
    view.cfl = Iteration::cfl;
    view.maxTimeRatio = nscom.max_time_ratio;
    view.cacheKey = grid;
    view.preferResidentState = true;
    view.downloadToHost = false;
    for ( int equation = 0; equation < view.nEquations; ++ equation )
        view.q[ equation ] = ( * unsf.q )[ equation ].data();
    HipFluxBackend::Shared().CalcLocalInviscidTimeStep( view );
#else
    throw std::runtime_error(
        "HIP local timestep was requested without HIP support" );
#endif
}

}

// Bind mean-flow (NS) pointers on global unsf for time-step kernels.
// Always pins NS_SOLVER (same contract as UNsField::Init).
void BindUNsFieldPointers()
{
    UnsGrid * grid = Zone::GetUnsGrid();
    unsf.invsr = ONEFLOW::GetFieldPointer< MRField >( grid, "invsr" );
    unsf.vissr = ONEFLOW::GetFieldPointer< MRField >( grid, "vissr" );
    unsf.gama  = ONEFLOW::GetFieldPointer< MRField >( grid, ZoneSharedFieldNames::Gamma );
    unsf.q  = FieldHome::GetSolutionField( grid, NS_SOLVER );
    unsf.q1 = FieldHome::GetStageSolutionField( grid, NS_SOLVER, 0 );
    unsf.q2 = FieldHome::GetStageSolutionField( grid, NS_SOLVER, 1 );
    // Update-delta database field; MULTI_STAGE may alias dq to residual in UNsCom::Init.
    unsf.dq = FieldHome::GetUpdateDeltaField( grid, NS_SOLVER );
    unsf.dqdx = FieldHome::GetGradField( grid, NS_SOLVER, 1 );
    unsf.dqdy = FieldHome::GetGradField( grid, NS_SOLVER, 2 );
    unsf.dqdz = FieldHome::GetGradField( grid, NS_SOLVER, 3 );

    // Shared zone physics fields (not solver-prefixed).
    unsf.visl = ONEFLOW::GetFieldPointer< MRField >( grid, ZoneSharedFieldNames::ViscosityLaminar );
    unsf.vist = ONEFLOW::GetFieldPointer< MRField >( grid, ZoneSharedFieldNames::ViscosityTurbulent );
    unsf.res = FieldHome::GetResidualField( grid, NS_SOLVER );
    unsf.res1 = FieldHome::GetStageResidualField( grid, NS_SOLVER, 0 );
    unsf.res2 = FieldHome::GetStageResidualField( grid, NS_SOLVER, 1 );
}

void InitTimeStepUns()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    BindUNsFieldPointers();

    FaceTopo & faceTopo = grid->GetFaceTopo();
    ug.lcf = & faceTopo.GetLeftCells();
    ug.rcf = & faceTopo.GetRightCells();

    FaceMesh & faceMesh = grid->GetFaceMesh();
    CellMesh & cellMesh = grid->GetCellMesh();

    ug.xfn = & faceMesh.xfn;
    ug.yfn = & faceMesh.yfn;
    ug.zfn = & faceMesh.zfn;
    ug.vfn = & faceMesh.vfn;
    ug.farea = & faceMesh.area;

    ug.vfx = & faceMesh.vfx;
    ug.vfy = & faceMesh.vfy;
    ug.vfz = & faceMesh.vfz;

    ug.xfc = & faceMesh.xfc;
    ug.yfc = & faceMesh.yfc;
    ug.zfc = & faceMesh.zfc;

    ug.xcc = & cellMesh.xcc;
    ug.ycc = & cellMesh.ycc;
    ug.zcc = & cellMesh.zcc;

    ug.cvol  = & cellMesh.vol;
    ug.cvol1 = & cellMesh.vol;
    ug.cvol2 = & cellMesh.vol;
    nscom.Init();
}

UTimeStep::UTimeStep()
{
    ;
}

UTimeStep::~UTimeStep()
{
    ;
}

void UTimeStep::Init()
{
    ug.Init();
    unsf.Init();
}

void UTimeStep::CalcTimeStep()
{
    this->Init();

    this->CalcCfl();

    const bool explicitRequest =
        HipProductionFeatureExplicitlyEnabled(
            "ONEFLOW_ENABLE_UNS_HIP_TIME_STEP" );
    bool useHipTimeStep = HipTimeStepRequested();
    if ( useHipTimeStep && ! HipTimeStepSupported() )
    {
        if ( explicitRequest )
        {
            throw std::runtime_error(
                "ONEFLOW_ENABLE_UNS_HIP_TIME_STEP=1 is unsupported for the "
                "current solver, timestep model, or accelerator chain" );
        }
        useHipTimeStep = false;
    }
    if ( useHipTimeStep )
    {
        CalcHipLocalInviscidTimeStep();
        return;
    }

    this->CalcSpectrumField();

    if ( nscom.timestepModel == 0 )
    {
        this->CalcLocalTimeStep();
    }
    else if ( nscom.timestepModel == 1 )
    {
        this->CalcGlobalTimeStep();
    }
    else
    {
        this->CalcLgTimeStep();
    }
}

void UTimeStep::CalcLocalTimeStep()
{
    this->CalcInvTimeStep();

    this->CalcVisTimeStep();

    this->ModifyTimeStep();
}

void UTimeStep::CalcInvTimeStep()
{
    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId  = cId;
        gcom.cvol  = ( * ug.cvol )[ cId ];
        nscom.invsr = ( * unsf.invsr )[ 0 ][ cId ];
        this->CalcCellInvTimeStep();
        ( * unsf.timestep )[ 0 ][ cId ] = nscom.timestep;
    }
}

void UTimeStep::CalcVisTimeStep()
{
    bool flag = vis_model.vismodel > 0 && nscom.visTimeStepModel > 0;
    if ( ! flag ) return;

    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ug.cId  = cId;
        gcom.cvol  = ( * ug.cvol )[ cId ];
        nscom.vissr = ( * unsf.vissr )[ 0 ][ cId ];
        nscom.timestep = ( * unsf.timestep )[ 0 ][ cId ];
        this->CalcCellVisTimeStep();
        ( * unsf.timestep )[ 0 ][ cId ] = nscom.timestep;
    }
}

void UTimeStep::CalcSpectrumField()
{
    this->CalcInvSpectrumField();
    this->CalcVisSpectrumField();
}

void UTimeStep::CalcInvSpectrumField()
{
    Grid & grid = Zone::GetGridReference();

    MRField * invsr = ONEFLOW::GetFieldPointer< MRField >( &grid, "invsr" );

    ONEFLOW::ZeroField( invsr, 1, grid.nCells );

    for ( int iFace = 0; iFace < grid.nFaces; ++ iFace )
    {
        this->SetId( iFace );

        this->PrepareData();

        this->CalcFaceInvSpec();

        this->UpdateInvSpectrumField();
    }
}

void UTimeStep::CalcVisSpectrumField()
{
    Grid & grid = Zone::GetGridReference();

    MRField * vissr = ONEFLOW::GetFieldPointer< MRField >( &grid, "vissr" );

    ONEFLOW::ZeroField( vissr, 1, grid.nCells );

    if ( vis_model.vismodel <= 0 ) return;

    for ( int iFace = 0; iFace < grid.nFaces; ++ iFace )
    {
        this->SetId( iFace );

        this->PrepareVisData();

        this->CalcFaceVisSpec();

        this->UpdateVisSpectrumField();
    }
}

void UTimeStep::SetId( int fId )
{
    ug.fId = fId;
    ug.lc = ( * ug.lcf )[ fId ];
    ug.rc = ( * ug.rcf )[ fId ];
}

void UTimeStep::PrepareData()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }

    nscom.gama1 = ( * unsf.gama  )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama  )[ 0 ][ ug.rc ];
    nscom.gama  = half * ( nscom.gama1 + nscom.gama2 );
}

void UTimeStep::PrepareVisData()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.rc ];
    }

    nscom.gama1 = ( * unsf.gama )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama )[ 0 ][ ug.rc ];

    gcom.cvol1 = ( * ug.cvol )[ ug.lc ];
    gcom.cvol2 = ( * ug.cvol )[ ug.rc ];

    nscom.visl1   = ( * unsf.visl  )[ 0 ][ ug.lc ];
    nscom.visl2   = ( * unsf.visl  )[ 0 ][ ug.rc ];

    nscom.visl = half * ( nscom.visl1 + nscom.visl2 );

    nscom.vist1 = ( * unsf.vist )[ 0 ][ ug.lc ];
    nscom.vist2 = ( * unsf.vist )[ 0 ][ ug.rc ];

    nscom.vist = half * ( nscom.vist1 + nscom.vist2 );

    gcom.xcc1 = ( * ug.xcc )[ ug.lc ];
    gcom.xcc2 = ( * ug.xcc )[ ug.rc ];

    gcom.ycc1 = ( * ug.ycc )[ ug.lc ];
    gcom.ycc2 = ( * ug.ycc )[ ug.rc ];

    gcom.zcc1 = ( * ug.zcc )[ ug.lc ];
    gcom.zcc2 = ( * ug.zcc )[ ug.rc ];
}

void UTimeStep::UpdateInvSpectrumField()
{
    ( * unsf.invsr )[ 0 ][ ug.lc ] += nscom.invsr;
    ( * unsf.invsr )[ 0 ][ ug.rc ] += nscom.invsr;
}

void UTimeStep::UpdateVisSpectrumField()
{
    ( * unsf.vissr )[ 0 ][ ug.lc ] += nscom.vissr;
    ( * unsf.vissr )[ 0 ][ ug.rc ] += nscom.vissr;
}

void UTimeStep::ModifyTimeStep()
{
    this->CalcMinTimeStep();
    if ( nscom.max_time_ratio <= 0 ) return;

    Real maxPermittedTimeStep = nscom.max_time_ratio * nscom.minTimeStep;
    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ( * unsf.timestep )[ 0 ][ cId ] = MIN( ( * unsf.timestep )[ 0 ][ cId ], maxPermittedTimeStep );
    }
}

void UTimeStep::CalcGlobalTimeStep()
{
    this->SetTimeStep( ctrl.pdt );
}

void UTimeStep::CalcMinTimeStep()
{
    nscom.minTimeStep = LARGE;
    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        nscom.minTimeStep = MIN( ( * unsf.timestep )[ 0 ][ cId ], nscom.minTimeStep );
    }
}

void UTimeStep::SetTimeStep( Real timestep )
{
    for ( int cId = 0; cId < ug.nCells; ++ cId )
    {
        ( * unsf.timestep )[ 0 ][ cId ] = timestep;
    }
}

void UTimeStep::CalcLgTimeStep()
{
    this->CalcLocalTimeStep();
    this->SetTimeStep( nscom.minTimeStep );
    ctrl.pdt = nscom.minTimeStep;
}


EndNameSpace
