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

#include "TurbInvFlux.h"
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
    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real vnl_p = half * ( vnl + ABS( vnl ) );
    Real vnr_n = half * ( vnr - ABS( vnr ) );

    Real rho_coef_l = inv.rho_coef * one + ( one - inv.rho_coef ) * inv.rl;
    Real rho_coef_r = inv.rho_coef * one + ( one - inv.rho_coef ) * inv.rr;

    Real coef_l = vnl_p * rho_coef_l;
    Real coef_r = vnr_n * rho_coef_r;

    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = coef_l * inv.prim1[ iEqu ] + coef_r * inv.prim2[ iEqu ];
    }
}


EndNameSpace