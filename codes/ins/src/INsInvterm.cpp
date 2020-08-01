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

//#include "UINsCorrectSpeed.h"
#include "INsInvterm.h"
#include "INsVisterm.h"
#include "UINsCom.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Ctrl.h"

BeginNameSpace( ONEFLOW )

INsInv iinv;



INsInv::INsInv()
{
    ;
}

INsInv::~INsInv()
{
    ;
}

void INsInv::Init()
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

INsInvterm::INsInvterm()
{
    ;
}

INsInvterm::~INsInvterm()
{
    ;
}

void INsInvterm::Solve()
{
}

void INsInvterm::CalcINsinvTerm()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	Real v2l = ONEFLOW::SQR(iinv.ul, iinv.vl, iinv.wl);
	Real v2r = ONEFLOW::SQR(iinv.ur, iinv.vr, iinv.wr);

	Real vnl = gcom.xfn * iinv.ul + gcom.yfn * iinv.vl + gcom.zfn * iinv.wl - gcom.vfn;       // V * n
	Real vnr = gcom.xfn * iinv.ur + gcom.yfn * iinv.vr + gcom.zfn * iinv.wr - gcom.vfn;

	Real rvnl = iinv.rl * vnl;   //ρ * V * n
	Real rvnr = iinv.rr * vnr;

	Real ratio = sqrt(iinv.rr / iinv.rl);
	Real coef = 1.0 / (1.0 + ratio);

	iinv.rm = (iinv.rl + iinv.rr ) * half;    //初始界面上的值（u、v、w ）
	iinv.um = (iinv.ul + iinv.ur ) * half;
	iinv.vm = (iinv.vl + iinv.vr ) * half;
	iinv.wm = (iinv.wl + iinv.wr ) * half;
	iinv.hm = (iinv.hl + iinv.hr * ratio) * coef;
	iinv.pm = (iinv.pl + iinv.pr * ratio) * coef;
	iinv.gama = (iinv.gama1 + iinv.gama2 * ratio) * coef;

	

	//Real v2 = ONEFLOW::SQR(iinv.um, iinv.vm, iinv.wm);


	iinv.vnflow = gcom.xfn * iinv.um + gcom.yfn * iinv.vm + gcom.zfn * iinv.wm;  //初始界面上 V*n

	iinv.fq0[ug.fId] = iinv.rl * iinv.vnflow * gcom.farea; //初始界面上的质量通量

	//Real gamm1 = iinv.gama - one;

	//Real c2 = gamm1 * (iinv.hm - half * v2);
	//iinv.cm = sqrt(ABS(c2));

	//iinv.aeig1 = ABS(iinv.vnrel);
	//iinv.aeig2 = ABS(iinv.vnrel + iinv.cm);
	//iinv.aeig3 = ABS(iinv.vnrel - iinv.cm);


	Real clr = MAX(0, iinv.fq0[ug.fId]);  //从界面左侧单元流入右侧单元的初始质量流量
	Real crl = clr - iinv.fq0[ug.fId];   //从界面右侧单元流入左侧单元的初始质量流量

	iinv.ai1[ug.lc] = iinv.ai1[ug.lc] + crl;   //界面左侧单元的系数
	iinv.ai2[ug.rc] = iinv.ai2[ug.rc] + clr;   //界面右侧单元的系数
	
	//iinv.flux[IIDX::IIRU] = iinv.rm * gcom.xfn * half* (iinv.ul + iinv.ur) * gcom.farea ;  
	//iinv.flux[IIDX::IIRV] = iinv.rm * gcom.yfn * half* (iinv.vl + iinv.vr)* gcom.farea ; 
	//iinv.flux[IIDX::IIRW] = iinv.rm * gcom.zfn * half* (iinv.vl + iinv.vr)* gcom.farea ; 
}

void INsInvterm::CalcINsFaceflux()
{
	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
	Real dy2 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
	Real dz2 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

	Real de1 = DIST(dx1, dy1, dz1);
	Real de2 = DIST(dx2, dy2, dz2);
	Real de = 1.0 / (de1 + de2);

	iinv.f1 = de2 * de;  //左单元权重
    iinv.f2 = de1 * de;  //右单元权重

	Real Va = iinv.f1[ug.fId]*(gcom.cvol1 / iinv.sp[ug.lc]) + iinv.f2[ug.fId] *(gcom.cvol2 / iinv.sp2[ug.rc]);   //（Vj/a）
    iinv.dist = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]); 
	Real Pd1 = visQ.dqdx1[IIDX::IIP] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + visQ.dqdy1[IIDX::IIP] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + visQ.dqdz1[IIDX::IIP] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);  //压力梯度项
	Real Pd2 = visQ.dqdx2[IIDX::IIP] * ((*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId]) + visQ.dqdy2[IIDX::IIP] * ((*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId]) + visQ.dqdz2[IIDX::IIP] * ((*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId]);
	Real Pd = Pd1 + Pd2;

	iinv.rm = (iinv.rl + iinv.rr) * half;  //界面密度
	iinv.um = (iinv.f1[ug.fId] *iinv.ul + iinv.f2[ug.fId] *iinv.ur) + (Va*gcom.xfn / iinv.dist[ug.fId])*( Pd - (iinv.pr - iinv.pl));  //界面密度
	iinv.vm = (iinv.f1[ug.fId] *iinv.vl + iinv.f2[ug.fId] *iinv.vr) + (Va*gcom.yfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));   //界面速度
	iinv.wm = (iinv.f1[ug.fId] *iinv.wl + iinv.f2[ug.fId] *iinv.wr) + (Va*gcom.zfn / iinv.dist[ug.fId])*(Pd - (iinv.pr - iinv.pl));

	iinv.flux[IIDX::IIRU] = iinv.rm * gcom.xfn * iinv.um * gcom.farea ;
	iinv.flux[IIDX::IIRV] = iinv.rm * gcom.yfn * iinv.vm * gcom.farea ;
	iinv.flux[IIDX::IIRW] = iinv.rm * gcom.zfn * iinv.wm * gcom.farea ;
}

void INsInvterm::CalcINsFaceCorrectPresscoef()
{
    iinv.Vdvj[ug.fId] = iinv.f1[ug.fId] * ( gcom.cvol1/((1+1)*iinv.sp[ug.lc]- iinv.spj[ug.lc])) + iinv.f2[ug.fId] * (gcom.cvol2 / ((1 + 1)*iinv.sp[ug.rc] - iinv.spj[ug.rc]));  // (Vp/dv)j，用于求面速度修正量
	iinv.aji[ug.fId] = iinv.rm * iinv.Vdvj[ug.fId] * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * gcom.farea / iinv.dist[ug.fId]; //ajp
}


  














EndNameSpace