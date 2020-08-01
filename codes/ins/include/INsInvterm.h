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

const int _ISCHEME_ROE            = 1;
const int _ISCHEME_VANLEER        = 2;
const int _ISCHEME_STEGER         = 3;
const int _ISCHEME_HLLE           = 4;
const int _ISCHEME_LAX_FRIEDRICHS = 5;
const int _ISCHEME_AUSMP          = 6;
const int _ISCHEME_AUSMPUP        = 7;
const int _ISCHEME_AUSMDV         = 8;
const int _ISCHEME_AUSMW          = 9;
const int _ISCHEME_AUSMPW         = 10;
const int _ISCHEME_HYBRIDROE      = 11;
const int _ISCHEME_SLAU2          = 12;


Real FMSplit1( const Real & mach, const Real & ipn );
Real FMSplit2( const Real & mach, const Real & ipn );
Real FMSplit4( const Real & mach, const Real & beta , const Real & ipn );
Real FPSplit5( const Real & mach, const Real & alpha, const Real & ipn );

class INsInv
{
public:
    INsInv();
    ~INsInv();
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
	RealField  ai1, ai2, bm, buc, bvc, bwc, bc,sp, sp1, sp2, spj, spu,spv,spw,fq0, aji, Vdv, Vdvj,bpu, bpv, bpw, spp, dist, f1, f2, ajp, app;
	RealField  pp, pp1,pp2,uu, vv, ww, uuj, vvj, wwj;//单元修正的压力和速度变量、界面修正速度变量
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

extern INsInv iinv;

class INsInvterm
{
public:
	INsInvterm();
    ~INsInvterm();
public:
    //typedef void (INsInvterm:: * InvtermPointer )();
public:
    void Solve();
public:
    //void SetPointer( int schemeIndex );
	//InvtermPointer InvtermPointer;
public:
	void CalcINsinvTerm();
	void CalcINsFaceflux();
	void CalcINsFaceCorrectPresscoef();
public:
};

EndNameSpace