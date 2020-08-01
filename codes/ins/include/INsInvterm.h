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
	//RealField * rf, *uf, *vf, *wf;
	RealField bc,aji, spj,Vdv, spp,  app, ajp,idx,spu,spv,spw,spu1, spv1, spw1, spu2, spv2, spw2, sp1, sp2, bm1,bm2, buc, bvc, bwc, aju, ajv, ajw, VdU, VdV, VdW,bpu, bpv, bpw, sppu, sppv, sppw, pp1, uu, vv, ww, ai1, ai2, aii1,aii2, akku1, akku2, akkv1, akkv2, akkw1, akkw2,aku1,aku2,akv1,akv2,akw1,akw2,f1, f2, pp, uuj, vvj, wwj,uc,vc,wc, mp,ppr,ppl, bp, dqqdx, dqqdy, dqqdz,muc,mvc,mwc,sju,sjv,sjw,pp0,pc,fq, up,vp,wp,spt, but, bvt,bwt,bmu1,bmu2,bmv1,bmv2,bmw1,bmw2, Fnu, Fnv, Fnw,  ukl, ukr, vkl, vkr, wkl,wkr, uml, umr, vml, vmr, wml, wmr, Pufd,Pvfd,Pwfd,Fpu,Fpv,Fpw,tux, tvy, twz, Fqu, Fqv, Fqw, FuT, FvT, FwT, PufdT, PvfdT, PwfdT, dsrl, elrn, visu,visv,visw, Fu1, Fv1, Fw1, bppu, bppv, bppw,ump,vmp,wmp, fqr,bi1,bi2,Fn, spc,fq1,fq2, Bpe1,Bpe2,mu,mv,mw,mpp, mua, mva, mwa,mppa, res_V, res_pp, res_up, res_vp, res_wp, rf, uf, vf, wf, Vdvu, Vdvv, Vdvw;
	RealField2D  spuj, spvj, spwj, sjp, ai, biu, biv,biw, sj,sd,sjd;
	

public:
    Real aeig1, aeig2, aeig3;
    Real meig1, meig2, meig3;

    Real eig11, eig12, eig13;
    Real eig21, eig22, eig23;
	Real vnrel, timestep, res_u, res_v, res_w, res_p, maxu, maxv, maxw, c2d, Fn1, Fn2, Fn3, lfdist, rfdist, ppd, vis, remax_V, remax_pp, remax_up, remax_vp, remax_wp, vnflow, l2rdx, l2rdy, l2rdz, Puf, Pvf, Pwf, Pdu, Pdv, Pdw, Ftu1, Ftv1, Ftw1, Ftu2, Ftv2, Ftw2, PufT, PvfT, PwfT, Pud, Pvd, Pwd, FtuT, FtvT, FtwT, Vau, Vav, Vaw, dist, Deun, Devn, Dewn, dlf, dfr, Bpe, fux, value;
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
	void CmpINsinvFlux();
	void CmpINsBcinvFlux();
	void CmpINsinvTerm();
	void CmpINsBcinvTerm();
	void CmpINsFaceflux();
	void CmpINsBcFaceflux();
	void CmpINsFaceCorrectPresscoef();
	void CmpINsBcFaceCorrectPresscoef();
	//void Roe      ();
    //void RoeOld   (){};
    //void HybridRoe(){};
    //void Vanleer  ();
    //void Steger   ();
    //void Hlle     ();
    //void LaxFriedrichs();
    //void Ausmp    ();
    //void AusmpUp  ();
    //void Ausmdv   ();
    //void Ausmw    ();
    //void Ausmpw   ();
    //void Slau2();
public:
   // void ModifyAbsoluteEigenvalue();
};

void INsCmpEnthalpy( RealField & prim, Real gama, Real & enthalpy );
void INsCmpTotalEnthalpyChange( RealField & prim, Real & gama, RealField & dq, Real & dh );
void INsPrimToQ( RealField & prim, Real gama, RealField & q );
void INsQToPrim( RealField & q, Real gama, RealField & prim, RealField & temp );
void INsCmpInternalEnergy( RealField & prim, Real gama, Real & em );

EndNameSpace