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

#include "ResidualBasePolicy.h"

#include "Ctrl.h"
#include "GridState.h"
#include "Multigrid.h"
#include "NsCom.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "TimeIntegral.h"
#include "HipProductionPolicy.h"

BeginNameSpace( ONEFLOW )

namespace
{

bool SolverSupportsDeviceRkBookkeeping()
{
    return HipProductionBackendSelected()
        && SolverState::solverType == NS_SOLVER
        && nscom.chemModel == 0
        && nscom.nTModel == 1
        && vis_model.vismodel == 0
        && ctrl.idualtime != 1;
}

}

bool DeviceOwnsExplicitRkBookkeeping()
{
#ifdef ONEFLOW_ENABLE_HIP
    const bool stateUpdate =
        HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UNS_HIP_STATE_UPDATE" );
    const bool reconstruction =
        HipProductionFeatureEnabled(
            "ONEFLOW_ENABLE_UNS_HIP_RECONSTRUCTION" );
    // The pure seven-argument truth table is unit tested; this call wires the
    // live globals into it so the test and production cannot drift apart.
    return SolverSupportsDeviceRkBookkeeping()
        && HipSkipHostRkSlotCopy(
               HipSkipHostRkSlotCopyRequested(),
               HipDeviceAuthoritativeStateEnabled(),
               stateUpdate && reconstruction,
               GridState::gridLevel == 0,
               MG::nMulti <= 1,
               ctrl.time_integral == MULTI_STAGE,
               ctrl.startStrategy != 2 );
#else
    return false;
#endif
}

bool HostResidualBookkeepingNeeded()
{
    // Coarse grids run the host update and read FIELD_FLOW/FIELD_RHS.
    if ( MG::nMulti > 1 ) return true;
    return ! DeviceOwnsExplicitRkBookkeeping();
}

EndNameSpace
