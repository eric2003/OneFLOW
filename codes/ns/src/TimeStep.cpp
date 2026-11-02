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

#include "TimeStep.h"
#include "TimeStepPolicy.h"
#include "NsInvFluxPolicy.h"
#include "Iteration.h"
#include "NsCom.h"
#include "HXMath.h"
#include "NsIdx.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

void TimeStep::CalcCfl()
{
    int iter = Iteration::outerSteps;

    if ( Iteration::dualtime == 1 )
    {
        iter = Iteration::innerSteps;
    }

    Iteration::cfl = Iteration::cfled;

    if ( iter < Iteration::ncfl )
    {
        Real delta = ( Iteration::cfled - Iteration::cflst ) / Iteration::ncfl;
        Iteration::cfl = Iteration::cflst + delta * iter;
    }
}

void TimeStep::CalcFaceInvSpec()
{
    Real rl = nscom.q1[ IDX::IR ];
    Real ul = nscom.q1[ IDX::IU ];
    Real vl = nscom.q1[ IDX::IV ];
    Real wl = nscom.q1[ IDX::IW ];
    Real pl = nscom.q1[ IDX::IP ];

    Real rr = nscom.q2[ IDX::IR ];
    Real ur = nscom.q2[ IDX::IU ];
    Real vr = nscom.q2[ IDX::IV ];
    Real wr = nscom.q2[ IDX::IW ];
    Real pr = nscom.q2[ IDX::IP ];

    Real vnl = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, ul, vl, wl );
    Real vnr = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, ur, vr, wr );

    nscom.gama = half * ( nscom.gama1 + nscom.gama2 );

    Real pm = half * ( pl + pr );
    Real rm = half * ( rl + rr );
    Real cm = IdealGasSoundSpeedFromPressure( rm, pm, nscom.gama );
    Real vn = half * ( vnl + vnr );

    nscom.invsr = InviscidFaceSpectralRadius( vn, cm, gcom.farea );
}

void TimeStep::CalcFaceVisSpec()
{
    const Real density = half * ( nscom.q1[ IDX::IR ] + nscom.q2[ IDX::IR ] );
    if ( nscom.visSRModel == 1 )
    {
        nscom.vissr = ViscousFaceSpectralRadiusModel1(
            density, nscom.visl, nscom.vist, nscom.oprl, nscom.oprt,
            nscom.gama, nscom.reynolds, gcom.farea );
    }
    else if ( nscom.visSRModel == 2 )
    {
        const Real dist = ABS(
            gcom.xfn * ( gcom.xcc2 - gcom.xcc1 )
            + gcom.yfn * ( gcom.ycc2 - gcom.ycc1 )
            + gcom.zfn * ( gcom.zcc2 - gcom.zcc1 ) );
        const Real viscosity = nscom.visl + nscom.vist;
        nscom.vissr = ViscousFaceSpectralRadiusModel2(
            density, viscosity, dist, nscom.reynolds, gcom.farea );
    }
    else if ( nscom.visSRModel == 3 )
    {
        nscom.vissr = ViscousFaceSpectralRadiusModel3(
            density, nscom.reynolds, gcom.farea );
    }
}

void TimeStep::CalcCellInvTimeStep()
{
    nscom.timestep = CellTimeStepFromSpectrum(
        Iteration::cfl, gcom.cvol, nscom.invsr );
}

void TimeStep::CalcCellVisTimeStep()
{
    const Real visTimeStep = CellTimeStepFromSpectrum(
        Iteration::cfl, gcom.cvol, nscom.vissr );
    nscom.timestep = CombineInvAndVisTimeStep( nscom.timestep, visTimeStep );
}

EndNameSpace
