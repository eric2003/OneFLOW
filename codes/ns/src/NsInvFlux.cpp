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

#include "NsInvFlux.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "HXMath.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

NsInv inv;

Real FMSplit1( const Real & mach, const Real & ipn )
{
    return half * ( mach + ipn * ABS( mach ) );
}

Real FMSplit2( const Real & mach, const Real & ipn )
{
    return ipn * fourth * SQR( mach + ipn );
}

Real FMSplit4( const Real & mach, const Real & beta, const Real & ipn )
{
    if ( ABS( mach ) >= one )
    {
        return FMSplit1( mach, ipn );
    }
    return FMSplit2( mach, ipn ) * ( one - ipn * sixteen * beta * FMSplit2( mach, - ipn ) );
}

Real FPSplit5( const Real & mach, const Real & alpha, const Real & ipn )
{
    Real macha = ABS( mach );
    if ( macha >= one )
    {
        return half * ( one + ipn * mach / ( macha + SMALL ) );
    }

    return FMSplit2( mach, ipn ) * ( ( two * ipn - mach ) - ipn * sixteen * alpha * mach * FMSplit2( mach,-ipn ) );
}

NsInv::NsInv()
{
    ;
}

NsInv::~NsInv()
{
    ;
}

void NsInv::Init()
{
    int nEqu = nscom.nEqu;
    prim.resize( nEqu );
    prim1.resize( nEqu );
    prim2.resize( nEqu );

    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    dq.resize( nEqu );

    flux.resize( nEqu );
    flux1.resize( nEqu );
    flux2.resize( nEqu );
}

NsInvFlux::NsInvFlux()
{
    ;
}

NsInvFlux::~NsInvFlux()
{
    ;
}

void NsInvFlux::SetPointer( int schemeId )
{
    invFluxPointer = & NsInvFlux::Roe;

    if ( schemeId == ISCHEME_ROE )
    {
        invFluxPointer = & NsInvFlux::Roe;
    }
    else if ( schemeId == ISCHEME_HYBRIDROE )
    {
        invFluxPointer = & NsInvFlux::HybridRoe;
    }
    else if ( schemeId == ISCHEME_AUSMP )
    {
        invFluxPointer = & NsInvFlux::Ausmp;
    }
    else if ( schemeId == ISCHEME_AUSMPUP )
    {
        invFluxPointer = & NsInvFlux::AusmpUp;
    }
    else if ( schemeId == ISCHEME_AUSMDV )
    {
        invFluxPointer = & NsInvFlux::Ausmdv;
    }
    else if ( schemeId == ISCHEME_AUSMW )
    {
        invFluxPointer = & NsInvFlux::Ausmw;
    }
    else if ( schemeId == ISCHEME_AUSMPW )
    {
        invFluxPointer = & NsInvFlux::Ausmpw;
    }
    else if ( schemeId == ISCHEME_VANLEER )
    {
        invFluxPointer = & NsInvFlux::Vanleer;
    }
    else if ( schemeId == ISCHEME_STEGER )
    {
        invFluxPointer = & NsInvFlux::Steger;
    }
    else if ( schemeId == ISCHEME_HLLE )
    {
        invFluxPointer = & NsInvFlux::Hlle;
    }
    else if ( schemeId == ISCHEME_LAX_FRIEDRICHS )
    {
        invFluxPointer = & NsInvFlux::LaxFriedrichs;
    }
    else if ( schemeId == ISCHEME_SLAU2 )
    {
        invFluxPointer = & NsInvFlux::Slau2;
    }
    else
    {
        invFluxPointer = & NsInvFlux::Roe;
    }
}

void NsInvFlux::Solve()
{
}

void NsInvFlux::ModifyAbsoluteEigenvalue()
{
    //Entropy fix
    inv.meig1 = inv.aeig1;
    inv.meig2 = inv.aeig2;
    inv.meig3 = inv.aeig3;

    if ( ctrl.ieigenfix == 1 )
    {
        Real maxEigenvalue = MAX( inv.aeig2, inv.aeig3 );
        Real m1 = maxEigenvalue * ctrl.centropy1;
        Real m2 = maxEigenvalue * ctrl.centropy2;

        inv.meig1 = MAX( m1, inv.aeig1 );
        inv.meig2 = MAX( m2, inv.aeig2 );
        inv.meig3 = MAX( m2, inv.aeig3 );
    }
}

void NsInvFlux::Roe()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real rvnl = inv.rl * vnl;
    Real rvnr = inv.rr * vnr;

    Real ratio = sqrt( inv.rr / inv.rl );
    Real coef  = 1.0 / ( 1.0 + ratio );

    inv.rm = sqrt( inv.rl * inv.rr );
    inv.um = ( inv.ul + inv.ur * ratio ) * coef;
    inv.vm = ( inv.vl + inv.vr * ratio ) * coef;
    inv.wm = ( inv.wl + inv.wr * ratio ) * coef;
    inv.hm = ( inv.hl + inv.hr * ratio ) * coef;
    inv.pm = ( inv.pl + inv.pr * ratio ) * coef;
    inv.gama = ( inv.gama1 + inv.gama2 * ratio ) * coef;

    Real v2 = ONEFLOW::SQR( inv.um, inv.vm, inv.wm );

    inv.vnflow = gcom.xfn * inv.um + gcom.yfn * inv.vm + gcom.zfn * inv.wm;
                 
    inv.vnrel = inv.vnflow - gcom.vfn;
                
    Real gamm1 = inv.gama - one;

    Real c2 = gamm1 * ( inv.hm - half * v2 );
    inv.cm = sqrt( ABS( c2 ) );

    inv.aeig1 = ABS( inv.vnrel                  );
    inv.aeig2 = ABS( inv.vnrel + inv.cm );
    inv.aeig3 = ABS( inv.vnrel - inv.cm );

    inv.flux1[ IDX::IR  ] = rvnl                            ;
    inv.flux1[ IDX::IRU ] = rvnl * inv.ul + gcom.xfn * inv.pl;
    inv.flux1[ IDX::IRV ] = rvnl * inv.vl + gcom.yfn * inv.pl;
    inv.flux1[ IDX::IRW ] = rvnl * inv.wl + gcom.zfn * inv.pl;
    inv.flux1[ IDX::IRE ] = rvnl * inv.hl + gcom.vfn * inv.pl;

    inv.flux2[ IDX::IR  ] = rvnr                            ;
    inv.flux2[ IDX::IRU ] = rvnr * inv.ur + gcom.xfn * inv.pr;
    inv.flux2[ IDX::IRV ] = rvnr * inv.vr + gcom.yfn * inv.pr;
    inv.flux2[ IDX::IRW ] = rvnr * inv.wr + gcom.zfn * inv.pr;
    inv.flux2[ IDX::IRE ] = rvnr * inv.hr + gcom.vfn * inv.pr;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux1[ iEqu ] = inv.prim1[ iEqu ] * inv.flux1[ IDX::IR ];
        inv.flux2[ iEqu ] = inv.prim2[ iEqu ] * inv.flux2[ IDX::IR ];
    }

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = inv.flux1[ iEqu ] + inv.flux2[ iEqu ];
    }

    PrimToQ( inv.prim1, inv.gama1, inv.q1 );
    PrimToQ( inv.prim2, inv.gama2, inv.q2 );

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.dq[ iEqu ] = inv.q2[ iEqu ] - inv.q1[ iEqu ];
    }

    this->ModifyAbsoluteEigenvalue();
               
    Real xi1 = ( two * inv.meig1 - inv.meig2 - inv.meig3 ) / ( two * c2 );
    Real xi2 = ( inv.meig2 - inv.meig3 ) / ( two * inv.cm );

    inv.prim[ IDX::IR ] = inv.rm;
    inv.prim[ IDX::IU ] = inv.um;
    inv.prim[ IDX::IV ] = inv.vm;
    inv.prim[ IDX::IW ] = inv.wm;
    inv.prim[ IDX::IP ] = inv.pm;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.prim[ iEqu ] = half * ( inv.prim1[ iEqu ] + inv.prim2[ iEqu ] );
    }

    Real dc   = inv.vnflow * inv.dq[ IDX::IR  ] - 
                gcom.xfn   * inv.dq[ IDX::IRU ] - 
                gcom.yfn   * inv.dq[ IDX::IRV ] - 
                gcom.zfn   * inv.dq[ IDX::IRW ];
    Real c2dc = c2 * dc;

    Real dh;
    CalcTotalEnthalpyChange( inv.prim, inv.gama, inv.dq, dh );

    Real term1 = dh   * xi1 + dc * xi2;
    Real term2 = c2dc * xi1 + dh * xi2;

    inv.flux[ IDX::IR  ] -= ( inv.meig1 * inv.dq[ IDX::IR  ] -          term1                      );
    inv.flux[ IDX::IRU ] -= ( inv.meig1 * inv.dq[ IDX::IRU ] - inv.um * term1 + gcom.xfn   * term2 );
    inv.flux[ IDX::IRV ] -= ( inv.meig1 * inv.dq[ IDX::IRV ] - inv.vm * term1 + gcom.yfn   * term2 );
    inv.flux[ IDX::IRW ] -= ( inv.meig1 * inv.dq[ IDX::IRW ] - inv.wm * term1 + gcom.zfn   * term2 );
    inv.flux[ IDX::IRE ] -= ( inv.meig1 * inv.dq[ IDX::IRE ] - inv.hm * term1 + inv.vnflow * term2 );

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] -= ( inv.meig1 * inv.dq[ iEqu ] - inv.prim[ iEqu ] * term1 );
    }

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] *= half;
    }
}

void NsInvFlux::Vanleer()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real rvnl = inv.rl * vnl;
    Real rvnr = inv.rr * vnr;

    Real c2l   = inv.gama1 * inv.pl / inv.rl;
    Real cl    = sqrt( ABS( c2l ) );
    Real machl = vnl / cl;

    Real c2r   = inv.gama2 * inv.pr / inv.rr;
    Real cr    = sqrt( ABS( c2r ) );
    Real machr = vnr / cr;

    if ( machl > one )
    {
        inv.flux1[ IDX::IR  ] = rvnl;
        inv.flux1[ IDX::IRU ] = rvnl * inv.ul + gcom.xfn * inv.pl;
        inv.flux1[ IDX::IRV ] = rvnl * inv.vl + gcom.yfn * inv.pl;
        inv.flux1[ IDX::IRW ] = rvnl * inv.wl + gcom.zfn * inv.pl;
        inv.flux1[ IDX::IRE ] = rvnl * inv.hl + gcom.vfn * inv.pl;
    }
    else if ( machl < - one )
    {
        inv.flux1[ IDX::IR  ] = zero;
        inv.flux1[ IDX::IRU ] = zero;
        inv.flux1[ IDX::IRV ] = zero;
        inv.flux1[ IDX::IRW ] = zero;
        inv.flux1[ IDX::IRE ] = zero;
    }
    else
    {
        Real fmass     = 0.25 * inv.rl * cl * SQR( machl + 1.0 );
        Real tmp       = ( - vnl + 2.0 * cl ) / inv.gama1;
        inv.flux1[ IDX::IR  ] = fmass;
        inv.flux1[ IDX::IRU ] = fmass * ( inv.ul + gcom.xfn * tmp );
        inv.flux1[ IDX::IRV ] = fmass * ( inv.vl + gcom.yfn * tmp );
        inv.flux1[ IDX::IRW ] = fmass * ( inv.wl + gcom.zfn * tmp );
        inv.flux1[ IDX::IRE ] = fmass * ( inv.hl + gcom.vfn * tmp );
    }

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux1[ iEqu ] = inv.prim1[ iEqu ] * inv.flux1[ IDX::IR ];
    }

    //Right contribution
    if ( machr < - one )
    {
        inv.flux2[ IDX::IR  ] = rvnr;
        inv.flux2[ IDX::IRU ] = rvnr * inv.ur + gcom.xfn * inv.pr;
        inv.flux2[ IDX::IRV ] = rvnr * inv.vr + gcom.yfn * inv.pr;
        inv.flux2[ IDX::IRW ] = rvnr * inv.wr + gcom.zfn * inv.pr;
        inv.flux2[ IDX::IRE ] = rvnr * inv.hr + gcom.vfn * inv.pr;
    }
    else if ( machr > one )
    {
        inv.flux2 = zero;
    }
    else
    {
        Real fmass     = - 0.25 * inv.rr * cr * SQR( machr - 1.0 );
        Real tmp       = ( - vnr - 2.0 * cr ) / inv.gama2;
        inv.flux2[ IDX::IR  ] = fmass;
        inv.flux2[ IDX::IRU ] = fmass * ( inv.ur + gcom.xfn * tmp );
        inv.flux2[ IDX::IRV ] = fmass * ( inv.vr + gcom.yfn * tmp );
        inv.flux2[ IDX::IRW ] = fmass * ( inv.wr + gcom.zfn * tmp );
        inv.flux2[ IDX::IRE ] = fmass * ( inv.hr + gcom.vfn * tmp );
    }

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux2[ iEqu ] = inv.prim2[ iEqu ] * inv.flux2[ IDX::IR ];
    }

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = inv.flux1[ iEqu ] + inv.flux2[ iEqu ];
    }
}

void NsInvFlux::Steger()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real rvnl = inv.rl * vnl;
    Real rvnr = inv.rr * vnr;

    inv.el  = inv.hl - inv.pl / ( inv.rl + SMALL );
    inv.er  = inv.hr - inv.pr / ( inv.rr + SMALL );

    Real c2l = inv.gama1 * inv.pl / inv.rl;
    Real c2r = inv.gama2 * inv.pr / inv.rr;

    Real cl  = sqrt( ABS( c2l ) );
    Real cr  = sqrt( ABS( c2r ) );

    inv.eig11 = vnl;
    inv.eig12 = vnl + cl;
    inv.eig13 = vnl - cl;

    inv.eig21 = vnr;
    inv.eig22 = vnr + cr;
    inv.eig23 = vnr - cr;

    Real aeig11 = ABS( inv.eig11 );
    Real aeig12 = ABS( inv.eig12 );
    Real aeig13 = ABS( inv.eig13 );

    Real aeig21 = ABS( inv.eig21 );
    Real aeig22 = ABS( inv.eig22 );
    Real aeig23 = ABS( inv.eig23 );

    Real mm = 0.01;

    Real m1 = MAX( ctrl.centropy1, mm );
    Real m2 = MAX( ctrl.centropy2, mm );

    //aeig11   = MAX( m1, aeig11 );
    //aeig12   = MAX( m2, aeig12 );
    //aeig13   = MAX( m2, aeig13 );

    //aeig21   = MAX( m1, aeig21 );
    //aeig22   = MAX( m2, aeig22 );
    //aeig23   = MAX( m2, aeig23 );

    aeig11   = DIST( aeig11, m1 );
    aeig12   = DIST( aeig12, m2 );
    aeig13   = DIST( aeig13, m2 );

    aeig21   = DIST( aeig21, m1 );
    aeig22   = DIST( aeig22, m2 );
    aeig23   = DIST( aeig23, m2 );

    inv.eig11 = half * ( inv.eig11 + aeig11 );
    inv.eig12 = half * ( inv.eig12 + aeig12 );
    inv.eig13 = half * ( inv.eig13 + aeig13 );

    inv.eig21 = half * ( inv.eig21 - aeig21 );
    inv.eig22 = half * ( inv.eig22 - aeig22 );
    inv.eig23 = half * ( inv.eig23 - aeig23 );

    Real c2gaml  = c2l / inv.gama1;
    Real c2gamr  = c2r / inv.gama2;

    Real term11 = c2gaml * ( two * inv.eig11 - inv.eig12 - inv.eig13 ) / ( c2l + c2l );
    Real term12 = c2gaml * ( inv.eig12 - inv.eig13 ) / ( cl + cl );

    Real term21 = c2gamr * ( two * inv.eig21 - inv.eig22 - inv.eig23 ) / ( c2r + c2r );
    Real term22 = c2gamr * ( inv.eig22 - inv.eig23 ) / ( cr + cr );

    Real vnFluid1 = vnl + gcom.vfn;
    Real vnFluid2 = vnr + gcom.vfn;

    inv.flux1[ IDX::IR  ] = inv.rl * ( inv.eig11          - term11                              );
    inv.flux1[ IDX::IRU ] = inv.rl * ( inv.eig11 * inv.ul - term11 * inv.ul + term12 * gcom.xfn );
    inv.flux1[ IDX::IRV ] = inv.rl * ( inv.eig11 * inv.vl - term11 * inv.vl + term12 * gcom.yfn );
    inv.flux1[ IDX::IRW ] = inv.rl * ( inv.eig11 * inv.wl - term11 * inv.wl + term12 * gcom.zfn );
    inv.flux1[ IDX::IRE ] = inv.rl * ( inv.eig11 * inv.el - term11 * inv.hl + term12 * vnFluid1 );

    inv.flux2[ IDX::IR  ] = inv.rr * ( inv.eig21          - term21                              );
    inv.flux2[ IDX::IRU ] = inv.rr * ( inv.eig21 * inv.ur - term21 * inv.ur + term22 * gcom.xfn );
    inv.flux2[ IDX::IRV ] = inv.rr * ( inv.eig21 * inv.vr - term21 * inv.vr + term22 * gcom.yfn );
    inv.flux2[ IDX::IRW ] = inv.rr * ( inv.eig21 * inv.wr - term21 * inv.wr + term22 * gcom.zfn );
    inv.flux2[ IDX::IRE ] = inv.rr * ( inv.eig21 * inv.er - term21 * inv.hr + term22 * vnFluid2 );

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux1[ iEqu ] = inv.prim1[ iEqu ] * inv.flux1[ IDX::IR ];
        inv.flux2[ iEqu ] = inv.prim2[ iEqu ] * inv.flux2[ IDX::IR ];
    }

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = inv.flux1[ iEqu ] + inv.flux2[ iEqu ];
    }
}

void NsInvFlux::Hlle()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real ratio = sqrt( inv.rr / inv.rl );
    Real coef  = 1.0 / ( 1.0 + ratio );

    inv.rm = sqrt( inv.rl * inv.rr );
    inv.um = ( inv.ul + inv.ur * ratio ) * coef;
    inv.vm = ( inv.vl + inv.vr * ratio ) * coef;
    inv.wm = ( inv.wl + inv.wr * ratio ) * coef;
    inv.hm = ( inv.hl + inv.hr * ratio ) * coef;
    inv.pm = ( inv.pl + inv.pr * ratio ) * coef;
    inv.gama = ( inv.gama1 + inv.gama2 * ratio ) * coef;
               
    Real gamm1 = inv.gama - one;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real rvnl = inv.rl * vnl;
    Real rvnr = inv.rr * vnr;

    inv.vnflow = gcom.xfn * inv.um + gcom.yfn * inv.vm + gcom.zfn * inv.wm;
    inv.vnrel = inv.vnflow - gcom.vfn;

    inv.cl  = sqrt( inv.gama1 * inv.pl / ( inv.rl + SMALL ) );
    inv.cr  = sqrt( inv.gama2 * inv.pr / ( inv.rr + SMALL ) );
    inv.cm  = sqrt( inv.gama * inv.pm / ( inv.rm + SMALL ) );

    inv.flux1[ IDX::IR  ] = rvnl                             ;
    inv.flux1[ IDX::IRU ] = rvnl * inv.ul + gcom.xfn * inv.pl;
    inv.flux1[ IDX::IRV ] = rvnl * inv.vl + gcom.yfn * inv.pl;
    inv.flux1[ IDX::IRW ] = rvnl * inv.wl + gcom.zfn * inv.pl;
    inv.flux1[ IDX::IRE ] = rvnl * inv.hl + gcom.vfn * inv.pl;

    inv.flux2[ IDX::IR  ] = rvnr                             ;
    inv.flux2[ IDX::IRU ] = rvnr * inv.ur + gcom.xfn * inv.pr;
    inv.flux2[ IDX::IRV ] = rvnr * inv.vr + gcom.yfn * inv.pr;
    inv.flux2[ IDX::IRW ] = rvnr * inv.wr + gcom.zfn * inv.pr;
    inv.flux2[ IDX::IRE ] = rvnr * inv.hr + gcom.vfn * inv.pr;

    PrimToQ( inv.prim1, inv.gama1, inv.q1 );
    PrimToQ( inv.prim2, inv.gama2, inv.q2 );

    Real bm = MIN( zero, MIN( inv.vnrel - inv.cm, vnl - inv.cl ) );
    Real bp = MAX( zero, MAX( inv.vnrel + inv.cm, vnr + inv.cr ) );

    Real c1 = bp * bm / ( bp - bm + SMALL );
    Real c2 = - half * ( bp + bm ) / ( bp - bm + SMALL );

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux1[ iEqu ] = inv.prim1[ iEqu ] * inv.flux1[ IDX::IR ];
        inv.flux2[ iEqu ] = inv.prim2[ iEqu ] * inv.flux2[ IDX::IR ];
    }

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = half * ( inv.flux1[ iEqu ] + inv.flux2[ iEqu ] )
                         + c1   * ( inv.q2   [ iEqu ] - inv.q1   [ iEqu ] )
                         + c2   * ( inv.flux2[ iEqu ] - inv.flux1[ iEqu ] );
    }
}

void NsInvFlux::LaxFriedrichs()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real ratio = sqrt( inv.rr / inv.rl );
    Real coef  = 1.0 / ( 1.0 + ratio );

    inv.rm = sqrt( inv.rl * inv.rr );
    inv.um = ( inv.ul + inv.ur * ratio ) * coef;
    inv.vm = ( inv.vl + inv.vr * ratio ) * coef;
    inv.wm = ( inv.wl + inv.wr * ratio ) * coef;
    inv.pm = ( inv.pl + inv.pr * ratio ) * coef;
    inv.gama = ( inv.gama1 + inv.gama2 * ratio ) * coef;

    Real gamm1 = inv.gama - one;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real rvnl = inv.rl * vnl;
    Real rvnr = inv.rr * vnr;

    inv.vnflow = gcom.xfn * inv.um + gcom.yfn * inv.vm + gcom.zfn * inv.wm;
    inv.vnrel = inv.vnflow - gcom.vfn;
    inv.cm  = sqrt( inv.gama * inv.pm / ( inv.rm + SMALL ) );

    inv.flux1[ IDX::IR  ] = rvnl                             ;
    inv.flux1[ IDX::IRU ] = rvnl * inv.ul + gcom.xfn * inv.pl;
    inv.flux1[ IDX::IRV ] = rvnl * inv.vl + gcom.yfn * inv.pl;
    inv.flux1[ IDX::IRW ] = rvnl * inv.wl + gcom.zfn * inv.pl;
    inv.flux1[ IDX::IRE ] = rvnl * inv.hl + gcom.vfn * inv.pl;

    inv.flux2[ IDX::IR  ] = rvnr                             ;
    inv.flux2[ IDX::IRU ] = rvnr * inv.ur + gcom.xfn * inv.pr;
    inv.flux2[ IDX::IRV ] = rvnr * inv.vr + gcom.yfn * inv.pr;
    inv.flux2[ IDX::IRW ] = rvnr * inv.wr + gcom.zfn * inv.pr;
    inv.flux2[ IDX::IRE ] = rvnr * inv.hr + gcom.vfn * inv.pr;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux1[ iEqu ] = inv.prim1[ iEqu ] * inv.flux1[ IDX::IR ];
        inv.flux2[ iEqu ] = inv.prim2[ iEqu ] * inv.flux2[ IDX::IR ];
    }

    PrimToQ( inv.prim1, inv.gama1, inv.q1 );
    PrimToQ( inv.prim2, inv.gama2, inv.q2 );

    Real maxeig = ABS( inv.vnrel ) + inv.cm;

    for ( int iEqu = 0; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = half * ( inv.flux1[ iEqu ] + inv.flux2[ iEqu ] ) - half * maxeig * ( inv.q2[ iEqu ] - inv.q1[ iEqu ] );
    }
}

void NsInvFlux::Ausmp()
{
    Real alphac = 3.0 / 16.0, betac = 0.125;

    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real orl = 1.0 / ( inv.rl + SMALL );
    Real orr = 1.0 / ( inv.rr + SMALL );

    inv.gama = half * ( inv.gama1 + inv.gama2 );

    Real c2l = inv.gama * inv.pl * orl;
    Real c2r = inv.gama * inv.pr * orr;

    inv.cl  = sqrt( c2l );
    inv.cr  = sqrt( c2r );
    inv.cm  = half * ( inv.cl + inv.cr );

    Real ocm = one / ( inv.cm + SMALL );

    Real ml =  vnl * ocm;
    Real mr =  vnr * ocm;

    Real fm4ml = FMSplit4( ml, betac,  one );
    Real fm4mr = FMSplit4( mr, betac, -one );

    Real fp5ml = FPSplit5( ml, alphac,  one );
    Real fp5mr = FPSplit5( mr, alphac, -one );

    Real mi  = fm4ml + fm4mr;
    Real p12 = fp5ml * inv.pl + fp5mr * inv.pr;

    Real mia =  ABS( mi );
    Real mpi =  half * ( mi + mia );
    Real mmi =  half * ( mi - mia );

    //AUSM + 
    inv.flux[ IDX::IR  ] = inv.cm * ( mpi * inv.rl          + mmi * inv.rr          );
    inv.flux[ IDX::IRU ] = inv.cm * ( mpi * inv.rl * inv.ul + mmi * inv.rr * inv.ur ) + gcom.xfn * p12;
    inv.flux[ IDX::IRV ] = inv.cm * ( mpi * inv.rl * inv.vl + mmi * inv.rr * inv.vr ) + gcom.yfn * p12;
    inv.flux[ IDX::IRW ] = inv.cm * ( mpi * inv.rl * inv.wl + mmi * inv.rr * inv.wr ) + gcom.zfn * p12;
    inv.flux[ IDX::IRE ] = inv.cm * ( mpi * inv.rl * inv.hl + mmi * inv.rr * inv.hr ) + gcom.vfn * p12;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = inv.cm * ( mpi * inv.rl * inv.prim1[ iEqu ] + mmi * inv.rr * inv.prim2[ iEqu ] );
    }
}

void NsInvFlux::AusmpUp()
{
    Real m2ref = SQR( nscom.mach_ref );

    Real fkp    = 0.25;
    Real fku    = 0.75;
    Real fsigma = one;

    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    inv.gama = half * ( inv.gama1 + inv.gama2 );

    Real orl = 1.0 / ( inv.rl + SMALL );
    Real orr = 1.0 / ( inv.rr + SMALL );

    Real c2l = inv.gama * inv.pl * orl;
    Real c2r = inv.gama * inv.pr * orr;

    //fourth speed interface
    inv.cl  = sqrt( c2l );
    inv.cr  = sqrt( c2r );
    inv.cm  = half * ( inv.cl + inv.cr );

    Real cm2 = inv.cm * inv.cm;
    Real ocm = one / ( inv.cm + SMALL );

    inv.rm  = half * ( inv.rl + inv.rr );

    Real ml  = vnl * ocm;
    Real mr  = vnr * ocm;

    Real mla = ABS( ml );
    Real mra = ABS( mr );

    Real ma2 = ( v2l + v2r ) / ( two * cm2 );
    Real m02 = MIN( one, MAX( ma2, m2ref ) );
    Real m0  = sqrt( m02 );

    Real fa  = m0 * ( two - m0 );

    Real alphac = 3.0 / 16.0 * ( - four + five * fa * fa );
    Real betac  = 0.125;

    Real fm4ml = FMSplit4( ml, betac,  one );
    Real fm4mr = FMSplit4( mr, betac, -one );

    Real fp5ml = FPSplit5( ml, alphac,  one );
    Real fp5mr = FPSplit5( mr, alphac, -one );

    Real mp  = - fkp / fa * MAX( one - fsigma * ma2, zero ) * ( inv.pr - inv.pl ) / ( inv.rm * cm2 );
    Real pu  = - fku * fp5ml * fp5mr * ( inv.rl + inv.rr ) * fa * inv.cm * ( vnr - vnl );

    //mp  = zero;
    pu  = zero;

    Real mi  = fm4ml + fm4mr + mp;
    Real p12 = fp5ml * inv.pl + fp5mr * inv.pr + pu;

    Real rvn;

    if ( mi > zero )
    {
        rvn = inv.cm * mi * inv.rl;
    }
    else
    {
        rvn = inv.cm * mi * inv.rr;
    }

    if ( rvn > zero )
    {
        inv.flux[ IDX::IR  ] = rvn                          ;
        inv.flux[ IDX::IRU ] = rvn * inv.ul + gcom.xfn * p12;
        inv.flux[ IDX::IRV ] = rvn * inv.vl + gcom.yfn * p12;
        inv.flux[ IDX::IRW ] = rvn * inv.wl + gcom.zfn * p12;
        inv.flux[ IDX::IRE ] = rvn * inv.hl + gcom.vfn * p12;

        for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
        {
            inv.flux[ iEqu ] = rvn * inv.prim1[ iEqu ];
        }
    }
    else
    {
        inv.flux[ IDX::IR  ] = rvn                          ;
        inv.flux[ IDX::IRU ] = rvn * inv.ur + gcom.xfn * p12;
        inv.flux[ IDX::IRV ] = rvn * inv.vr + gcom.yfn * p12;
        inv.flux[ IDX::IRW ] = rvn * inv.wr + gcom.zfn * p12;
        inv.flux[ IDX::IRE ] = rvn * inv.hr + gcom.vfn * p12;

        for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
        {
            inv.flux[ iEqu ] = rvn * inv.prim2[ iEqu ];
        }
    }
}

void NsInvFlux::Ausmdv()
{
    Real alphac = 3.0 / 16.0, betac = 0.125;
    Real ssw = 0.0, ssw_a;
    Real K_SWITCH = 10.0, C_EFIX = 0.125;

    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    inv.gama = half * ( inv.gama1 + inv.gama2 );

    Real orl = 1.0 / ( inv.rl + SMALL );
    Real orr = 1.0 / ( inv.rr + SMALL );

    Real c2l = inv.gama * inv.pl * orl;
    Real c2r = inv.gama * inv.pr * orr;

    //这个公式应该参照文献好好对一下，很可能在运动网格下有问题！！！！！！！eric 20140126
    //这里的vn_l和vn_r是惯性坐标系下的法向速度
    //而vnl和vnr是相对于运动界面的法向速度
    Real vn_l  = vnl + gcom.vfn;
    Real vn_r  = vnr + gcom.vfn;

    //fourth speed interface
    inv.cl  = sqrt( c2l );
    inv.cr  = sqrt( c2r );
    inv.cm  = MAX( inv.cl, inv.cr );
    Real ocm = one / ( inv.cm + SMALL );

    ssw = zero;
    //Determine shock switch
    if ( ( vnl > inv.cl && vnr < inv.cr ) || ( vnl > - inv.cl && vnr < - inv.cr ) )
    {
        ssw = one;
    }
    ssw_a = one - ssw;

    //Haenel / van Leer
    Real ml    =  vnl / inv.cl;
    Real mla   =  ABS( ml );

    Real mr    =  vnr / inv.cr;
    Real mra   =  ABS( mr );

    Real mp1 = half * ( ml + mla );
    Real mm1 = half * ( mr - mra );

    Real mpl, ppl;

    if ( mla >= one )
    {
        mpl = mp1;
        ppl = half * ( one + ml / ( mla + SMALL ) );
    }
    else
    {
        Real tmp = fourth * SQR( ml + one );

        mpl = tmp;
        ppl = tmp * ( two - ml );
    }

    Real mmr, pmr;

    if ( mra >= one )
    {
        mmr = mm1;
        pmr = half * ( one - mr / ( mra + SMALL ) );
    }
    else
    {
        Real tmp = fourth * SQR( mr - one );

        mmr = - tmp;
        pmr = tmp * ( two + mr );
    }

    Real p12_v =  ppl * inv.pl + pmr * inv.pr;

    Real rvnl_v = inv.cl * inv.rl * mpl;
    Real rvnr_v = inv.cr * inv.rr * mmr;

    //AUSMDV
    ml    =  vnl * ocm;
    mla   =  ABS( ml );

    mr    =  vnr * ocm;
    mra   =  ABS( mr );

    mp1 = half * ( ml + mla );
    mm1 = half * ( mr - mra );

    Real r1  = inv.pl * orl;
    Real r2  = inv.pr * orr;
    Real r3  = two / ( r1 + r2 );
    Real alphal = r1 * r3;
    Real alphar = r2 * r3;

    if ( mla >= one )
    {
        mpl = mp1;
        ppl = half * ( one + ml / ( mla + SMALL ) );
    }
    else
    {
        Real tmp = fourth * SQR( ml + one );

        mpl = alphal * tmp + ( one - alphal ) * mp1;
        ppl = tmp * ( two - ml );
    }

    if ( mra >= one )
    {
        mmr = mm1;
        pmr = half * ( one - mr / ( mra + SMALL ) );
    }
    else
    {
        Real tmp = fourth * SQR( mr - one );

        mmr = - alphar * tmp + ( one - alphar ) * mm1;
        pmr = tmp * ( two + mr );
    }

    Real p12_d =  ppl * inv.pl + pmr * inv.pr;

    Real rvnl_d = inv.cm * inv.rl * mpl;
    Real rvnr_d = inv.cm * inv.rr * mmr;

    Real rvn    = rvnl_d + rvnr_d;
    rvnl_d = half * ( rvn + ABS( rvn ) );
    rvnr_d = half * ( rvn - ABS( rvn ) );

    //rvnl = ssw * rvnl_v + ssw_a * rvnl_d;
    //rvnr = ssw * rvnr_v + ssw_a * rvnr_d;

    Real rvnl = rvnl_v;
    Real rvnr = rvnr_v;
    Real p12  = p12_v ;

    //AUSMD/AUSMV switch
    Real s = half * MIN( one, K_SWITCH * ABS( inv.pr - inv.pl ) / MIN( inv.pl, inv.pr ) );

    Real rvvn = ssw * ( rvnl_v * vn_l + rvnr_v * vn_r );
    rvvn += ( half - s ) * ssw_a *          (   rvnl_d * vn_l +   rvnr_d * vn_r )
          + ( half + s ) * ssw_a * inv.cm * ( inv.rl * mpl * vn_l + inv.rr * mmr * vn_r );

    //entropy fix
    bool atmp = ( ( vnl - inv.cl ) < 0.0 ) && ( ( vnr - inv.cr ) > 0.0 );
    bool btmp = ( ( vnl + inv.cl ) < 0.0 ) && ( ( vnr + inv.cr ) > 0.0 );

    //atmp = btmp = false;

    if ( ( atmp && ( ! btmp ) ) || ( ( ! atmp ) && btmp ) )
    {
        if ( atmp && ( ! btmp ) )
        {
            Real tmp   = ssw_a * C_EFIX * ( ( vnr - inv.cr ) - ( vnl - inv.cl ) );
            rvnl += tmp * inv.rl;
            rvnr -= tmp * inv.rr;
        }
        else
        {
            Real tmp   = ssw_a * C_EFIX * ( ( vnr + inv.cr ) - ( vnl + inv.cl ) );
            rvnl += tmp * inv.rl;
            rvnr -= tmp * inv.rr;
        }
    }

    Real dul = inv.ul - gcom.xfn * vn_l;
    Real dvl = inv.vl - gcom.yfn * vn_l;
    Real dwl = inv.wl - gcom.zfn * vn_l;

    Real dur = inv.ur - gcom.xfn * vn_r;
    Real dvr = inv.vr - gcom.yfn * vn_r;
    Real dwr = inv.wr - gcom.zfn * vn_r;

    inv.flux[ IDX::IR  ] = rvnl           + rvnr      ;
    inv.flux[ IDX::IRU ] = rvnl * dul     + rvnr * dur;
    inv.flux[ IDX::IRV ] = rvnl * dvl     + rvnr * dvr;
    inv.flux[ IDX::IRW ] = rvnl * dwl     + rvnr * dwr;
    inv.flux[ IDX::IRE ] = rvnl * inv.hl  + rvnr * inv.hr;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = rvnl * inv.prim1[ iEqu ] + rvnr * inv.prim2[ iEqu ];
    }

    rvvn += p12;

    inv.flux[ IDX::IRU ] += rvvn * gcom.xfn;
    inv.flux[ IDX::IRV ] += rvvn * gcom.yfn;
    inv.flux[ IDX::IRW ] += rvvn * gcom.zfn;
    inv.flux[ IDX::IRE ] += p12  * gcom.vfn;
}

void NsInvFlux::Ausmw()
{
    Real alphac = 3.0 / 16.0, betac = 0.125;

    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = ONEFLOW::SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = ONEFLOW::SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    inv.gama = half * ( inv.gama1 + inv.gama2 );

    Real orl = 1.0 / ( inv.rl + SMALL );
    Real orr = 1.0 / ( inv.rr + SMALL );

    Real vnla  =  ABS( vnl );
    Real alc   =  inv.gama * inv.pl * orl;

    Real vnra  =  ABS( vnr );
    Real arc   =  inv.gama * inv.pr * orr;

    Real al = sqrt( alc );
    Real ar = sqrt( arc );
    Real ai = half * ( al + ar );

    Real rai =  one / ( ai + SMALL );

    Real ml  =  vnl * rai;
    Real mla =  ABS( ml );

    Real mr  =  vnr * rai;
    Real mra =  ABS( mr );

    Real mpl = half * ( ml + mla );
    Real ppl = half * ( one + ml / ( mla + SMALL ) );

    if ( mla < one )
    {
        Real tmp1 = ml + one;
        Real tmp2 = tmp1 * tmp1;
        Real tmp3 = ( ml * ml - one );
        Real tmp4 = tmp3 * tmp3;

        mpl = fourth * tmp2 + betac * tmp4;
        ppl = fourth * tmp2 * ( two - ml ) + alphac * ml * tmp4;
    }

    Real mml = ml  - mpl;
    Real pml = one - ppl;

    Real mpr = half * ( mr + mra );
    Real ppr = half * ( one + mr / ( mra + SMALL ) );

    if ( mra < one )
    {
        Real tmp1 = mr + one;
        Real tmp2 = tmp1 * tmp1;
        Real tmp3 = ( mr * mr - one );
        Real tmp4 = tmp3 * tmp3;

        mpr = fourth * tmp2 + betac * tmp4;
        ppr = fourth * tmp2 * ( two - mr ) + alphac * mr * tmp4;
    }

    Real mmr = mr - mpr;
    Real pmr = one - ppr;

    Real presi  = ppl * inv.pl + pmr * inv.pr;

    //!此mi相当于m1/2，与文献定义一致
    Real mi  = mpl + mmr;

    Real mia = ABS( mi );

    Real mpi = half * ( mi + mia );
    Real mmi = half * ( mi - mia );

    //AUSM + _W
    Real rtr = one / ( inv.rl + inv.rr + SMALL );
    Real sp  = two * inv.rr * rtr;
    Real sm  = two * inv.rl * rtr;

    Real ccl = sp * mpl + ( one - sp ) * mpi;
    Real ccr = sm * mmr + ( one - sm ) * mmi;

    inv.flux[ IDX::IR  ] = ai * ( mpi * inv.rl          + mmi * inv.rr          );
    inv.flux[ IDX::IRU ] = ai * ( ccl * inv.rl * inv.ul + ccr * inv.rr * inv.ur ) + gcom.xfn * presi;
    inv.flux[ IDX::IRV ] = ai * ( ccl * inv.rl * inv.vl + ccr * inv.rr * inv.vr ) + gcom.yfn * presi;
    inv.flux[ IDX::IRW ] = ai * ( ccl * inv.rl * inv.wl + ccr * inv.rr * inv.wr ) + gcom.zfn * presi;
    inv.flux[ IDX::IRE ] = ai * ( mpi * inv.rl * inv.hl + mmi * inv.rr * inv.hr ) + gcom.vfn * presi;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = ai * ( mpi * inv.rl * inv.prim1[ iEqu ] + mmi * inv.rr * inv.prim2[ iEqu ] );
    }
}

void NsInvFlux::Ausmpw()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real orl = 1.0 / ( inv.rl + SMALL );
    Real orr = 1.0 / ( inv.rr + SMALL );

    inv.gama = half * ( inv.gama1 + inv.gama2 );

    Real c2l = inv.gama1 * inv.pl * orl;
    Real c2r = inv.gama2 * inv.pr * orr;
    inv.cl  = sqrt( c2l );
    inv.cr  = sqrt( c2r );
    inv.cm  = half * ( inv.cl + inv.cr );

    Real ocm = one / ( inv.cm + SMALL );

    Real ml  =  vnl * ocm;
    Real mr  =  vnr * ocm;

    Real fm4ml = FMSplit4( ml, zero,  one );
    Real fm4mr = FMSplit4( mr, zero, -one );

    Real fp5ml = FPSplit5( ml, zero,  one );
    Real fp5mr = FPSplit5( mr, zero, -one );

    Real m12  = fm4ml + fm4mr;
    Real p12  = fp5ml * inv.pl + fp5mr * inv.pr;

    Real minplr = MIN( inv.pl / ( inv.pr + SMALL ), inv.pr / ( inv.pl + SMALL ) );
    Real fw     = one - minplr * minplr * minplr;

    Real op12 = one / ( p12 + SMALL );

    Real fl = zero;
    Real fr = zero;

    if ( ABS( ml ) < one )
    {
        fl = inv.pl * op12 - one;
    }

    if ( ABS( mr ) < one )
    {
        fr = inv.pr * op12 - one;
    }

    Real mpl, mmr;

    if ( m12 >= zero )
    {
        Real oneFr = one + fr;
        mpl = fm4ml + fm4mr * ( ( one - fw ) * oneFr - fl );
        mmr = fm4mr * fw * oneFr;
    }
    else
    {
        Real oneFl = one + fl;
        mpl = fm4ml * fw * oneFl;
        mmr = fm4mr + fm4ml * ( ( one - fw ) * oneFl - fr );
    }

    inv.flux[ IDX::IR  ] = inv.cm * ( mpl * inv.rl          + mmr * inv.rr          );
    inv.flux[ IDX::IRU ] = inv.cm * ( mpl * inv.rl * inv.ul + mmr * inv.rr * inv.ur ) + gcom.xfn * p12;
    inv.flux[ IDX::IRV ] = inv.cm * ( mpl * inv.rl * inv.vl + mmr * inv.rr * inv.vr ) + gcom.yfn * p12;
    inv.flux[ IDX::IRW ] = inv.cm * ( mpl * inv.rl * inv.wl + mmr * inv.rr * inv.wr ) + gcom.zfn * p12;
    inv.flux[ IDX::IRE ] = inv.cm * ( mpl * inv.rl * inv.hl + mmr * inv.rr * inv.hr ) + gcom.vfn * p12;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = inv.cm * ( mpl * inv.rl * inv.prim1[ iEqu ] + mmr * inv.rr * inv.prim2[ iEqu ] );
    }
}

void NsInvFlux::Slau2()
{
    Extract( inv.prim1, inv.rl, inv.ul, inv.vl, inv.wl, inv.pl );
    Extract( inv.prim2, inv.rr, inv.ur, inv.vr, inv.wr, inv.pr );

    Real v2l = SQR( inv.ul, inv.vl, inv.wl );
    Real v2r = SQR( inv.ur, inv.vr, inv.wr );

    Real hint_l, hint_r;

    CalcEnthalpy( inv.prim1, inv.gama1, hint_l );
    CalcEnthalpy( inv.prim2, inv.gama2, hint_r );

    inv.hl  = hint_l + half * v2l;
    inv.hr  = hint_r + half * v2r;

    Real vnl  = gcom.xfn * inv.ul + gcom.yfn * inv.vl + gcom.zfn * inv.wl - gcom.vfn;
    Real vnr  = gcom.xfn * inv.ur + gcom.yfn * inv.vr + gcom.zfn * inv.wr - gcom.vfn;

    Real orl = 1.0 / ( inv.rl + SMALL );
    Real orr = 1.0 / ( inv.rr + SMALL );

    Real c2l = inv.gama1 * inv.pl * orl;
    Real c2r = inv.gama2 * inv.pr * orr;

    inv.cl  = sqrt( c2l );
    inv.cr  = sqrt( c2r );

    inv.cm  = half * ( inv.cl + inv.cr ); //middle sound speed    

    Real t = inv.rr / inv.rl;
    Real vna = ( ABS( vnl ) + t * ABS( vnr ) ) / ( 1 + t );

    Real ml = vnl / inv.cm;
    Real mr = vnr / inv.cm;

    Real g = - MAX( MIN( ml, 0.0), -1.0 ) * MIN( MAX( mr, 0.0), 1.0 );

    Real vnp = ( 1 - g ) * vna + g * ABS( vnl );
    Real vnn = ( 1 - g ) * vna + g * ABS( vnr );

    Real va = sqrt( half * ( v2l + v2r ) );
    Real m12 = MIN( 1.0, va / inv.cm );
    Real ka = SQR( 1 - m12 );

    Real ms = half * ( inv.rl * ( vnl + vnp ) + inv.rr * ( vnr - vnn ) - ka * ( inv.pr - inv.pl ) / inv.cm );

    Real fp5ml = FPSplit5( ml, zero,  one );
    Real fp5mr = FPSplit5( mr, zero, -one );

    Real ps = 0;
    ps += half * ( inv.pl + inv.pr );
    ps += half * ( fp5ml - fp5mr ) * ( inv.pl - inv.pr );
    ps += va *( fp5ml + fp5mr - 1 ) * sqrt( inv.rl * inv.rr ) * inv.cm;

    Real mp = half * ( ms + ABS( ms ) );
    Real mn = half * ( ms - ABS( ms ) );

    inv.flux[ IDX::IR  ] = ( mp          + mn          );
    inv.flux[ IDX::IRU ] = ( mp * inv.ul + mn * inv.ur ) + gcom.xfn * ps;
    inv.flux[ IDX::IRV ] = ( mp * inv.vl + mn * inv.vr ) + gcom.yfn * ps;
    inv.flux[ IDX::IRW ] = ( mp * inv.wl + mn * inv.wr ) + gcom.zfn * ps;
    inv.flux[ IDX::IRE ] = ( mp * inv.hl + mn * inv.hr ) + gcom.vfn * ps;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        inv.flux[ iEqu ] = ( mp * inv.prim1[ iEqu ] + mn * inv.prim2[ iEqu ] );
    }
}

void CalcEnthalpy( RealField & prim, Real gama, Real & enthalpy )
{
    Real & density  = prim[ IDX::IR ];
    Real & pressure = prim[ IDX::IP ];

    if ( nscom.chemModel == 0 )
    {
        enthalpy = ( gama / ( gama - one ) ) * ( pressure / density );
    }
    else
    {
    }
}

void CalcTotalEnthalpyChange( RealField & prim, Real & gama, RealField & dq, Real & dh )
{
    Real v2, ae, af;

    Real & density  = prim[ IDX::IR ];
    Real & um       = prim[ IDX::IU ];
    Real & vm       = prim[ IDX::IV ];
    Real & wm       = prim[ IDX::IW ];
    Real & pressure = prim[ IDX::IP ];

    v2 = ONEFLOW::SQR( um, vm, wm );

    ae = gama - one;
    af = half * ae * v2;
    dh = - ae * ( um * dq[ IDX::IRU ] + vm * dq[ IDX::IRV ] + wm * dq[ IDX::IRW ] - dq[ IDX::IRE ] );
 
    dh += af * dq[ IDX::IR  ];
}

void PrimToQ( RealField & prim, Real gama, RealField & q )
{
    Real em;
    Real & density  = prim[ IDX::IR ];
    Real &      um  = prim[ IDX::IU ];
    Real &      vm  = prim[ IDX::IV ];
    Real &      wm  = prim[ IDX::IW ];
    Real & pressure = prim[ IDX::IP ];
    
    ONEFLOW::CalcInternalEnergy( prim, gama, em );

    q[ IDX::IR  ] = density;
    q[ IDX::IRU ] = density * um;
    q[ IDX::IRV ] = density * vm;
    q[ IDX::IRW ] = density * wm;
    q[ IDX::IRE ] = density * em;

    for ( int iEqu = nscom.nBEqu; iEqu < nscom.nEqu; ++ iEqu )
    {
        q[ iEqu ] = density * prim[ iEqu ];
    }
}

void QToPrim( RealField & q, Real gama, RealField & prim, RealField & temp )
{
    Real density = ABS( q[ IDX::IR ] );
    prim[ IDX::IR ] = density;

    Real rd = 1.0 / density;
    Real um  = q[ IDX::IU ] * rd;
    Real vm  = q[ IDX::IV ] * rd;
    Real wm  = q[ IDX::IW ] * rd;
    Real rem = q[ IDX::IP ];
    Real v2  = ONEFLOW::SQR( um, vm, wm );

    Real reint = rem - half * density * v2;

    prim[ IDX::IR ] = density;
    prim[ IDX::IU ] = um;
    prim[ IDX::IV ] = vm;
    prim[ IDX::IW ] = wm;

    Real pressure = ( gama - one ) * reint;
    prim[ IDX::IP ] = pressure;
}


void CalcInternalEnergy( RealField & prim, Real gama, Real & em )
{
    Real & rm = prim[ IDX::IR ];
    Real & um = prim[ IDX::IU ];
    Real & vm = prim[ IDX::IV ];
    Real & wm = prim[ IDX::IW ];
    Real & pm = prim[ IDX::IP ];
    Real v2 = ONEFLOW::SQR( um, vm, wm );

    em = ( pm / rm ) / ( gama - one ) + half * v2;
}


EndNameSpace