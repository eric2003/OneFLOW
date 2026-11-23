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

#pragma once
#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

// One decision shared by the two places that must agree about the explicit
// Runge-Kutta bookkeeping:
//   * FieldTaskReg decides whether LOAD_Q / LOAD_RESIDUALS still fill the
//     host FIELD_FLOW / FIELD_RHS slots;
//   * UNsInvFlux decides whether the device residual starts from the uploaded
//     host field or from zero on the device.
// If only one of them changes, the stage silently accumulates on the previous
// stage's residual.  That is a numerical regression, not a performance knob.
//
// The device owns this bookkeeping only for the validated explicit chain:
// steady multi-stage, finest grid, single grid level, NS solver, inviscid,
// no chemistry, no dual time.  Everything else keeps the host slots
// authoritative (multigrid restriction reads FIELD_FLOW, the host update and
// the multigrid defect read FIELD_RHS).
//
// Implementation: codes/solver/src/ResidualBasePolicy.cpp.
bool DeviceOwnsExplicitRkBookkeeping();

// The named residual field and the FIELD_RHS slot are the host-side residual
// bookkeeping used by the multigrid defect (restrict/prolong reads FIELD_RHS)
// and by the host explicit update.  When the device owns both the base and the
// accumulation for a single-grid steady step, nothing reads them: the device
// residual is the only one that matters and it is downloaded only for output.
bool HostResidualBookkeepingNeeded();

EndNameSpace
