/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "ITimestep.h"
#include "Iteration.h"
#include "INsCom.h"
#include "HXMath.h"
#include "INsIDX.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

ITimestep::ITimestep()
{
    ;
}

ITimestep::~ITimestep()
{
    ;
}

void ITimestep::CmpCfl()
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

void ITimestep::CmpFaceInvSpec()
{
    Real rl = inscom.q1[ IIDX::IIR ];
    Real ul = inscom.q1[ IIDX::IIU ];
    Real vl = inscom.q1[ IIDX::IIV ];
    Real wl = inscom.q1[ IIDX::IIW ];
    Real pl = inscom.q1[ IIDX::IIP ];
    Real cl = sqrt( inscom.gama1 * pl / rl );

    Real rr = inscom.q2[ IIDX::IIR ];
    Real ur = inscom.q2[ IIDX::IIU ];
    Real vr = inscom.q2[ IIDX::IIV ];
    Real wr = inscom.q2[ IIDX::IIW ];
    Real pr = inscom.q2[ IIDX::IIP ];
    Real cr = sqrt( inscom.gama2 * pr / rr );
        
    Real vnl  = gcom.xfn * ul + gcom.yfn * vl + gcom.zfn * wl - gcom.vfn;
    Real vnr  = gcom.xfn * ur + gcom.yfn * vr + gcom.zfn * wr - gcom.vfn;

    inscom.gama = half * ( inscom.gama1 + inscom.gama2 );

    Real pm = half * ( pl + pr );
    Real rm = half * ( rl + rr );
    Real cm = sqrt( inscom.gama * pm / rm );
    Real vn = half * ( vnl + vnr );

    inscom.invsr = half * gcom.farea * ( ABS( vn ) + cm );
}

void ITimestep::CmpFaceVisSpec()
{
    if ( inscom.visSRModel == 1 )
    {
        Real density = half * ( inscom.q1[ IIDX::IIR ] + inscom.q2[ IIDX::IIR ] );

        Real c1 = 4.0 / 3.0 * ( inscom.visl + inscom.vist );
        Real c2 = inscom.gama * ( inscom.visl * inscom.oprl + inscom.vist * inscom.oprt );
        Real c3 = two * MAX( c1, c2 ) / ( inscom.reynolds * density );
        Real farea2 = SQR( gcom.farea );

        inscom.vissr = farea2 * c3;
    }
    else if ( inscom.visSRModel == 2 )
    {
        Real dist = ABS(  gcom.xfn * ( gcom.xcc2 - gcom.xcc1 )
                        + gcom.yfn * ( gcom.ycc2 - gcom.ycc1 )
                        + gcom.zfn * ( gcom.zcc2 - gcom.zcc1 ) );

        Real viscosity = inscom.visl + inscom.vist;
        Real density   = half * ( inscom.q1[ IIDX::IIR ] + inscom.q2[ IIDX::IIR ] );

        Real c1  = 2.0 * viscosity / ( density * dist * inscom.reynolds + SMALL );
        inscom.vissr = half * c1 * gcom.farea;
    }
    else if ( inscom.visSRModel == 3 )
    {
        Real density = half * ( inscom.q1[ IIDX::IIR ] + inscom.q2[ IIDX::IIR ] );
        Real farea2 = SQR( gcom.farea );

        inscom.vissr = farea2 / ( inscom.reynolds * density + SMALL );
    }
}

void ITimestep::CmpCellInvTimestep()
{
    inscom.timestep = Iteration::cfl * gcom.cvol / inscom.invsr;
}

void ITimestep::CmpCellVisTimestep()
{
    Real visTimestep = Iteration::cfl * gcom.cvol / inscom.vissr;
    inscom.timestep *= visTimestep / ( inscom.timestep + visTimestep );
}

EndNameSpace