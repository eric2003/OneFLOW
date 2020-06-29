/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "Iteration.h"
#include "NsCom.h"
#include "HXMath.h"
#include "NsIdx.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

TimeStep::TimeStep()
{
    ;
}

TimeStep::~TimeStep()
{
    ;
}

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
    Real cl = sqrt( nscom.gama1 * pl / rl );

    Real rr = nscom.q2[ IDX::IR ];
    Real ur = nscom.q2[ IDX::IU ];
    Real vr = nscom.q2[ IDX::IV ];
    Real wr = nscom.q2[ IDX::IW ];
    Real pr = nscom.q2[ IDX::IP ];
    Real cr = sqrt( nscom.gama2 * pr / rr );
        
    Real vnl  = gcom.xfn * ul + gcom.yfn * vl + gcom.zfn * wl - gcom.vfn;
    Real vnr  = gcom.xfn * ur + gcom.yfn * vr + gcom.zfn * wr - gcom.vfn;

    nscom.gama = half * ( nscom.gama1 + nscom.gama2 );

    Real pm = half * ( pl + pr );
    Real rm = half * ( rl + rr );
    Real cm = sqrt( nscom.gama * pm / rm );
    Real vn = half * ( vnl + vnr );

    nscom.invsr = half * gcom.farea * ( ABS( vn ) + cm );
}

void TimeStep::CalcFaceVisSpec()
{
    if ( nscom.visSRModel == 1 )
    {
        Real density = half * ( nscom.q1[ IDX::IR ] + nscom.q2[ IDX::IR ] );

        Real c1 = 4.0 / 3.0 * ( nscom.visl + nscom.vist );
        Real c2 = nscom.gama * ( nscom.visl * nscom.oprl + nscom.vist * nscom.oprt );
        Real c3 = two * MAX( c1, c2 ) / ( nscom.reynolds * density );
        Real farea2 = SQR( gcom.farea );

        nscom.vissr = farea2 * c3;
    }
    else if ( nscom.visSRModel == 2 )
    {
        Real dist = ABS(  gcom.xfn * ( gcom.xcc2 - gcom.xcc1 )
                        + gcom.yfn * ( gcom.ycc2 - gcom.ycc1 )
                        + gcom.zfn * ( gcom.zcc2 - gcom.zcc1 ) );

        Real viscosity = nscom.visl + nscom.vist;
        Real density   = half * ( nscom.q1[ IDX::IR ] + nscom.q2[ IDX::IR ] );

        Real c1  = 2.0 * viscosity / ( density * dist * nscom.reynolds + SMALL );
        nscom.vissr = half * c1 * gcom.farea;
    }
    else if ( nscom.visSRModel == 3 )
    {
        Real density = half * ( nscom.q1[ IDX::IR ] + nscom.q2[ IDX::IR ] );
        Real farea2 = SQR( gcom.farea );

        nscom.vissr = farea2 / ( nscom.reynolds * density + SMALL );
    }
}

void TimeStep::CalcCellInvTimeStep()
{
    nscom.timestep = Iteration::cfl * gcom.cvol / nscom.invsr;
}

void TimeStep::CalcCellVisTimeStep()
{
    Real visTimeStep = Iteration::cfl * gcom.cvol / nscom.vissr;
    nscom.timestep *= visTimeStep / ( nscom.timestep + visTimeStep );
}

EndNameSpace