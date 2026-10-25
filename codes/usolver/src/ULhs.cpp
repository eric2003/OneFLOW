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

#include "ULhs.h"
#include "Zone.h"
#include "Ctrl.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "UNsCom.h"
#include "UCom.h"
#include "FieldWrap.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "NsCom.h"
#include "AccelRuntime.h"
#include "AccelViews.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipFluxBackend.h"
#endif
#include <cstdlib>
#include <stdexcept>
#include <iostream>


BeginNameSpace( ONEFLOW )

namespace
{

bool UseHipDeviceStateUpdate( int solverType )
{
    const char * feature = "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE";
    const bool explicitRequest =
        HipProductionFeatureExplicitlyEnabled( feature );
    if ( ! HipProductionFeatureEnabled( feature ) )
        return false;
    if ( solverType != NS_SOLVER )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update is only supported for the NS solver" );
    }
    if ( ! HipProductionFeatureEnabled(
             "ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION" ) )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update requires HIP device reconstruction" );
    }
    if ( nscom.chemModel != 0 || nscom.nTModel != 1 )
    {
        if ( ! explicitRequest ) return false;
        throw std::runtime_error(
            "HIP state update requires a single ideal-gas temperature field" );
    }
#ifndef ONEFLOW_ENABLE_HIP
    throw std::runtime_error(
        "HIP state update was requested without HIP support" );
#else
    return true;
#endif
}

bool UseResidentHipTimeStep()
{
    // Follow the backend policy: the device timestep is not mirrored back to
    // the host, so CALC_LHS must consume the resident one whenever the
    // resident timestep stage is active.
    return HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_TIME_STEP" );
}

}

void ULhs::CalcLHS( int solverType )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    if ( UseHipDeviceStateUpdate( solverType ) )
    {
#ifdef ONEFLOW_ENABLE_HIP
        ug.Init();
        unsf.Init();
        CellStateUpdateView view;
        view.nCells = ug.nCells;
        view.nGhostCells = ug.nTCell - ug.nCells;
        view.nEquations = nscom.nEqu;
        view.timeStep = ( * unsf.timestep )[ 0 ].data();
        view.cellVolume = ug.cvol->data();
        view.gamma = nscom.gama_ref;
        view.rkCoefficient = ctrl.lhscoef;
        view.cacheKey = grid;
        view.preferResidentTimeStep = UseResidentHipTimeStep();
        for ( int equation = 0; equation < view.nEquations; ++ equation )
        {
            view.primitive[ equation ] = ( * unsf.q )[ equation ].data();
        }
        HipFluxBackend::Shared().ScaleCurrentResidual( view );
        return;
#endif
    }

    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    // LHS scales the residual for the given solverType (NS "res" / Turb "turbres"),
    // not the database update-delta field. Under MULTI_STAGE, *Com may alias dq
    // to this same residual storage.
    MRField * residual = FieldHome::GetResidualField( grid, solverType );

    for ( int iEqu = 0; iEqu < solverInfo->nTEqu; ++ iEqu )
    {
        for ( int cId = 0; cId < ug.nCells; ++ cId )
        {
            //Real dt = ctrl.pdt;
            Real dt = ( * unsf.timestep )[ 0 ][ cId ];
            ( * residual )[ iEqu ][ cId ] *= ( dt / ( * ug.cvol )[ cId ] ) * ctrl.lhscoef;
        }
    }
}

EndNameSpace
