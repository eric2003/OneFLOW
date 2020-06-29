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


#pragma once
#include "HXDefine.h"
#include "HXArray.h"
BeginNameSpace( ONEFLOW )

const int ISCHEME_ROE            = 1;
const int ISCHEME_VANLEER        = 2;
const int ISCHEME_STEGER         = 3;
const int ISCHEME_HLLE           = 4;
const int ISCHEME_LAX_FRIEDRICHS = 5;
const int ISCHEME_AUSMP          = 6;
const int ISCHEME_AUSMPUP        = 7;
const int ISCHEME_AUSMDV         = 8;
const int ISCHEME_AUSMW          = 9;
const int ISCHEME_AUSMPW         = 10;
const int ISCHEME_HYBRIDROE      = 11;
const int ISCHEME_SLAU2          = 12;


Real FMSplit1( const Real & mach, const Real & ipn );
Real FMSplit2( const Real & mach, const Real & ipn );
Real FMSplit4( const Real & mach, const Real & beta , const Real & ipn );
Real FPSplit5( const Real & mach, const Real & alpha, const Real & ipn );

class NsInv
{
public:
    NsInv();
    ~NsInv();
public:
    void Init();
public:
    Real gama;
    Real gama1;
    Real gama2;
public:
    RealField prim, prim1, prim2;
    RealField q, q1, q2;
    RealField dq;
    RealField flux, flux1, flux2;
public:
    Real aeig1, aeig2, aeig3;
    Real meig1, meig2, meig3;

    Real eig11, eig12, eig13;
    Real eig21, eig22, eig23;

    Real vnrel, vnflow;
    Real cl, cr, cm;
public:
    Real rl, ul, vl, wl, pl, hl, el;
    Real rr, ur, vr, wr, pr, hr, er;
    Real rm, um, vm, wm, pm, hm, em;
};

extern NsInv inv;

class NsInvFlux
{
public:
    NsInvFlux ();
    ~NsInvFlux();
public:
    typedef void ( NsInvFlux:: * InvFluxPointer )();
public:
    void Solve();
public:
    void SetPointer( int schemeIndex );
    InvFluxPointer invFluxPointer;
public:
    void Roe      ();
    void RoeOld   (){};
    void HybridRoe(){};
    void Vanleer  ();
    void Steger   ();
    void Hlle     ();
    void LaxFriedrichs();
    void Ausmp    ();
    void AusmpUp  ();
    void Ausmdv   ();
    void Ausmw    ();
    void Ausmpw   ();
    void Slau2();
public:
    void ModifyAbsoluteEigenvalue();
};

void CalcEnthalpy( RealField & prim, Real gama, Real & enthalpy );
void CalcTotalEnthalpyChange( RealField & prim, Real & gama, RealField & dq, Real & dh );
void PrimToQ( RealField & prim, Real gama, RealField & q );
void QToPrim( RealField & q, Real gama, RealField & prim, RealField & temp );
void CalcInternalEnergy( RealField & prim, Real gama, Real & em );

EndNameSpace