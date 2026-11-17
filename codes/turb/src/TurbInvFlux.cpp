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

#include "TurbInvFlux.h"
#include "TurbInvFluxPolicy.h"
#include "NsInvFluxPolicy.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "HXMath.h"
#include "NsCtrl.h"

BeginNameSpace( ONEFLOW )

TurbInv turbInv;

TurbInv::TurbInv()
{
    ;
}

TurbInv::~TurbInv()
{
    ;
}

void TurbInv::Init()
{
    int nEqu = 1;
    if ( vis_model.vismodel == 3 )
    {
        nEqu = 1;
        rho_coef = 1;
    }
    else if ( vis_model.vismodel == 4 )
    {
        nEqu = 2;
        rho_coef = 0;
    }

    prim1.resize( nEqu );
    prim2.resize( nEqu );

    q1.resize( nEqu );
    q2.resize( nEqu );

    flux.resize( nEqu );
}

TurbInvFlux::TurbInvFlux()
{
    ;
}

TurbInvFlux::~TurbInvFlux()
{
    ;
}

void TurbInvFlux::RoeFlux()
{
    TurbInv & inv = turbInv;
    int nEqu = inv.prim1.size();
    const Real vnl = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, inv.ul, inv.vl, inv.wl );
    const Real vnr = RelativeFaceNormalVelocity(
        gcom.xfn, gcom.yfn, gcom.zfn, gcom.vfn, inv.ur, inv.vr, inv.wr );

    Real vnl_p, vnl_n, vnr_p, vnr_n;
    SplitRelativeNormalVelocity( vnl, vnl_p, vnl_n );
    SplitRelativeNormalVelocity( vnr, vnr_p, vnr_n );
    (void) vnl_n;
    (void) vnr_p;

    const Real rho_coef_l = TurbConvectiveDensityWeight( inv.rho_coef, inv.rl );
    const Real rho_coef_r = TurbConvectiveDensityWeight( inv.rho_coef, inv.rr );

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = UpwindScalarConvectiveFlux(
            vnl_p, rho_coef_l, inv.prim1[ iEqu ],
            vnr_n, rho_coef_r, inv.prim2[ iEqu ] );
    }
}


EndNameSpace
