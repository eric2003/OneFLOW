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
#include "TurbTrans.h"

BeginNameSpace( ONEFLOW )

const int ISA = 0;
const int IKE = 0;
const int IKW = 1;
const int ITGama   = 2;
const int ITRect   = 3;
const int ITReGama = 2;

const int DES   = 1;
const int DDES  = 2;
const int IDDES = 3;

class TurbCom
{
public:
    TurbCom();
    ~TurbCom();
public:
    bool init_flag;
    int nEqu;
    Real visl1, visl2;
    Real vist1, vist2;
    Real visl, vist, vis;
    Real rho, rho1, rho2;

    Real reynolds, oreynolds;
    Real trans_df, trans_dct;
    Transition trans;
public:
    int sst_type;
    int ibld;
    int transition_model;
    int des_model;
    int iprod_sa, iprod_sst;
    int ft2_flag;
    int iturb_visflux;
    int easm_model;
    Real rprod;
    Real max_vis_ratio;
    Real ref_sa, ref_sst;
    Real inflow_intensity, inflow_viscosity;
    Real turb_cfl_ratio;
public:
    Real dudx, dudy, dudz;
    Real dvdx, dvdy, dvdz;
    Real dwdx, dwdy, dwdz;
    Real work1, work2, work3, work4;
    Real str, stp;
    Real compress;
    Real nuet;
    Real dqdxSa, dqdySa, dqdzSa;
    Real len_scale, dist;
    Real spec_sa, res_sa;
public:
    Real cdes, cdes_ke, cdes_kw;
public:
    Real xsi, xsi2, xsi3;
    Real sac2, sac3;
    Real karm, karm2, okarm2;
    Real a1;
    Real cb1, cb2;
    Real cv1, cv13;
    Real cw1, cw2, cw3, cw36;
    Real fv1, fv2;
    Real xk, xw, xk2, fk, fw, fwStar;
    Real ct3, ct4;
    Real cb2s, cw1k;
    Real ft2, ftrans;
public:
    Real pklim, fbetas, fbeta, betas;
    Real cmu, clim;
public:
    Real bld, bld1, bld2;
    Real sigk, sigk1, sigk2;
    Real sigw, sigw1, sigw2;
    Real sigd, sigd1, sigd2;
    Real alphaw, alphaw1, alphaw2;
    Real beta, beta1, beta2;

    Real srck, srcw;
    Real diak, diaw;
    Real fskn, fswn;
    Real oork, oorw;
    Real ke, kw;
    Real dkedx, dkedy, dkedz;
    Real dkwdx, dkwdy, dkwdz;
    Real kelim, kwlim;
    Real kGamaLim, kRectLim;
    Real crelax;
    int nneg, npos, nmax;
    Real maxvist;
    Real maxvistratio;

    int maxid;

    Real srcg, srcr;
    Real speg, sper;
public:
    Real um, vm, wm;
    Real intermittency, rectabar;
public:
    Real s11, s22, s33;
    Real s12, s13, s23;
    Real w12, w13, w23;
    Real gmeff, sij2, divv;
    Real prodk, dissk;
    Real prodw, dissw, cdkww;
    Real cross_term;
    Real crossdiff, cdkwmin, cdkwmax;
public:
    Real sigma, osigma;
    RealField comVis;
    RealField flux;
public:
    int bctype;
    int bcdtkey;
    int isowallbc;
    int turb_ilim, tns_ilim;
public:
    int faceOuterNormal;

    RealField prims1;
    RealField prims2;

    RealField primt1;
    RealField primt2;

    RealField inflow;
    RealField * bcflow;

    RealField ns_prims1;
    RealField ns_prims2;

    RealField q, q0, dq;
public:
    void Init();
    void InitConst();
    void InitInflow();
    void CalcSigkw();
    void CalcWorkVar();
    void CalcSaProd();
    void CalcVGrad();
    void CalcProdk();
    void CalcDissk();
    void LimitProdk();
    void CalcProdwKwMenter();
    void CalcProdwKwWilcox1998();
    void CalcProdwKwWilcox2006();
    void CalcProdwKwDefault();
    void CalcProdwEasmKw2003();
    void CalcProdwEasmKw2005();
    void ModifyPd();
    void CalcSrc();
    void CalcCrossDiff();
    void CalcCdkwmin();
    void CalcCellBlendingTerm();
    void CalcCrossing();
    void CalcFbetaOfKwWilcox1998();
    void CalcFbetaOfKwWilcox2006();
    void CalcFbetaOfEasmKw2003();
    void RGamaTransition();
    void CalcSrcSa();
};

extern TurbCom turbcom;

EndNameSpace