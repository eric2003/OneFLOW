/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "Iteration.h"
#include "UINsCom.h"
#include "Zone.h"
#include "DataBase.h"
#include "UCom.h"
#include "Com.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Ctrl.h"
#include "Boundary.h"
#include "BcRecord.h"

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

void INsInvterm::CalcINsinvFlux()
{

		INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

		INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);



		iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //Values at the initial interface (u, v, w)

		iinv.uf[ug.fId] = (iinv.ul + iinv.ur) * half;

		iinv.vf[ug.fId] = (iinv.vl + iinv.vr) * half;

		iinv.wf[ug.fId] = (iinv.wl + iinv.wr) * half;

		iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;  //V * n at initial interface

		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //Mass flux at the initial interface

}

void INsInvterm::CalcINsBcinvFlux()
{

	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);

	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	int bcType = ug.bcRecord->bcType[ug.fId];


	iinv.rf[ug.fId] = (iinv.rl + iinv.rr) * half;    //Values at the initial interface (u, v, w)

	iinv.uf[ug.fId] = (iinv.ul + iinv.ur) * half;

	iinv.vf[ug.fId] = (iinv.vl + iinv.vr) * half;

	iinv.wf[ug.fId] = (iinv.wl + iinv.wr) * half;

	iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId] - gcom.vfn;  //V * n at initial interface

	if(bcType == BC::SOLID_SURFACE)
	{
		iinv.fq[ug.fId] = 0;
	}
	else
	{
		iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * gcom.farea; //Mass flux at the initial interface
	}
}

void INsInvterm::CalcINsinvTerm()
{ 
		Real clr = MAX(0, iinv.fq[ug.fId]);  //Mass flow from the left cell of the interface to the right cell

		Real crl = clr - iinv.fq[ug.fId];   //Mass flow from the right cell of the interface to the left cell
		
		iinv.ai[ug.fId][0] = clr;
		iinv.ai[ug.fId][1] = crl;

		/*iinv.ai[0][ug.fId] = clr;
		iinv.ai[1][ug.fId] = crl;*/

}

void INsInvterm::CalcINsBcinvTerm()
{
	
		Real clr = MAX(0, iinv.fq[ug.fId]);  //Initial mass flow rate from the left cell of the interface to the right cell
		Real crl = clr - iinv.fq[ug.fId];   //Initial mass flow rate from the right cell of the interface to the left cell


		iinv.aii1[ug.fId] = crl;   //The flow rate of the left cell of the surface flow
		iinv.aii2[ug.fId] = clr;   //The flow rate of the right cell of the surface flow

		iinv.ai1[ug.lc] += crl;   //Flow into the unit
		iinv.ai2[ug.rc] += clr;   //Flow out of the unit
}

void INsInvterm::CalcINsFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	
	iinv.Vau = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.xfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.xfn / (iinv.spc[ug.rc])); //Df*n
	iinv.Vav = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.yfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.yfn / (iinv.spc[ug.rc]));
	iinv.Vaw = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] * gcom.zfn / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] * gcom.zfn / (iinv.spc[ug.rc]));

	iinv.dist = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.Deun = iinv.Vau / iinv.dist;   //Df*n/e*n
	iinv.Devn = iinv.Vav / iinv.dist;
	iinv.Dewn = iinv.Vaw / iinv.dist;


	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xcc)[ug.rc] - (*ug.xfc)[ug.fId];
	Real dy2 = (*ug.ycc)[ug.rc] - (*ug.yfc)[ug.fId];
	Real dz2 = (*ug.zcc)[ug.rc] - (*ug.zfc)[ug.fId];

	iinv.dlf = sqrt(dx1*dx1 + dy1 * dy1 + dz1 * dz1);
	iinv.dfr = sqrt(dx2*dx2 + dy2 * dy2 + dz2 * dz2);

	iinv.Bpe = (*uinsf.dqdx)[IIDX::IIP][ug.lc] * (dx1) + (*uinsf.dqdy)[IIDX::IIP][ug.lc] * (dy1) + (*uinsf.dqdz)[IIDX::IIP][ug.lc] * (dz1) +
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] * (dx2) + (*uinsf.dqdy)[IIDX::IIP][ug.rc] * (dy2) + (*uinsf.dqdz)[IIDX::IIP][ug.rc] * (dz2) -
		(iinv.pr - iinv.pl);

	iinv.rf[ug.fId] = half * (iinv.rl+ iinv.rr);
	iinv.uf[ug.fId] = (iinv.f1[ug.fId] * iinv.ul + iinv.f2[ug.fId] * iinv.ur)+iinv.Deun * iinv.Bpe;  //Interface prediction speed at the next moment
	iinv.vf[ug.fId] = (iinv.f1[ug.fId] * iinv.vl + iinv.f2[ug.fId] * iinv.vr)+iinv.Devn * iinv.Bpe;
	iinv.wf[ug.fId] = (iinv.f1[ug.fId] * iinv.wl + iinv.f2[ug.fId] * iinv.wr)+iinv.Dewn * iinv.Bpe;
	

	iinv.vnflow = (*ug.xfn)[ug.fId] * iinv.uf[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf[ug.fId];

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * (*ug.farea)[ug.fId];  //Next moment interface prediction flux


	Real clr = MAX(0, iinv.fq[ug.fId]);  //Initial mass flow rate from the left cell of the interface to the right cell

	Real crl = clr - iinv.fq[ug.fId];   //Initial mass flow rate from the right cell of the interface to the left cell

	iinv.ai[ug.fId][0] = clr;
	iinv.ai[ug.fId][1] = crl;

}


void INsInvterm::CalcINsBcFaceflux()
{
	INsExtract(iinv.prim1, iinv.rl, iinv.ul, iinv.vl, iinv.wl, iinv.pl);
	INsExtract(iinv.prim2, iinv.rr, iinv.ur, iinv.vr, iinv.wr, iinv.pr);

	iinv.dist = gcom.xfn * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + gcom.yfn * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + gcom.zfn * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);
	
	int bcType = ug.bcRecord->bcType[ug.fId];

	iinv.rf[ug.fId] = (iinv.rl + iinv.rr)*half;
	iinv.uf[ug.fId] = (iinv.ul + iinv.ur)*half;
	iinv.vf[ug.fId] = (iinv.vl + iinv.vr)*half;
	iinv.wf[ug.fId] = (iinv.wl + iinv.wr)*half;
	
	iinv.vnflow = gcom.xfn * iinv.uf[ug.fId] + gcom.yfn * iinv.vf[ug.fId] + gcom.zfn * iinv.wf[ug.fId];

	iinv.fq[ug.fId] = iinv.rf[ug.fId] * iinv.vnflow * (*ug.farea)[ug.fId];  //Next moment interface prediction flux
	

	/*iinv.uf1[ug.fId] = (iinv.ul + iinv.ur)*half;
	iinv.vf1[ug.fId] = (iinv.vl + iinv.vr)*half;
	iinv.wf1[ug.fId] = (iinv.wl + iinv.wr)*half;

	iinv.uf2[ug.fId] = 0;
	iinv.vf2[ug.fId] = 0;
	iinv.wf2[ug.fId] = 0;

	iinv.uf[ug.fId] = iinv.uf1[ug.fId] + iinv.uf2[ug.fId];
	iinv.vf[ug.fId] = iinv.vf1[ug.fId] + iinv.vf2[ug.fId];
	iinv.wf[ug.fId] = iinv.wf1[ug.fId] + iinv.wf2[ug.fId];

	iinv.fq1[ug.fId] = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uf1[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf1[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf1[ug.fId])*(*ug.farea)[ug.fId];
	iinv.fq2[ug.fId] = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uf2[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vf2[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wf2[ug.fId])*(*ug.farea)[ug.fId];
	iinv.fq[ug.fId] = iinv.fq1[ug.fId] + iinv.fq2[ug.fId];*/


	if (bcType == BC::SOLID_SURFACE)
	{
		iinv.fq[ug.fId] = 0;
	}

	Real clr = MAX(0, iinv.fq[ug.fId]);  //Initial mass flow rate from the left cell of the interface to the right cell

	Real crl = clr - iinv.fq[ug.fId];   //Initial mass flow rate from the right cell of the interface to the left cell

	iinv.ai[ug.fId][0] = clr;
	iinv.ai[ug.fId][1] = crl;
									  
}

void INsInvterm::CalcINsFaceCorrectPresscoef()
{
	iinv.Vdvu[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] /(iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]  / (iinv.spc[ug.rc]));  // -MF * n, used to calculate the surface velocity correction
	iinv.Vdvv[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc] / (iinv.spc[ug.rc]));
	iinv.Vdvw[ug.fId] = iinv.f1[ug.fId] * ((*ug.cvol1)[ug.lc] / (iinv.spc[ug.lc])) + iinv.f2[ug.fId] * ((*ug.cvol2)[ug.rc]  / (iinv.spc[ug.rc]));
	

	/*if ((*ug.xfn)[ug.fId] >= 0)
	{
		iinv.Vdvu[ug.fId] =  ((*ug.farea)[ug.fId] / (iinv.spc[ug.lc]));
	}
	else
	{
		iinv.Vdvu[ug.fId] = ((*ug.farea)[ug.fId] / (iinv.spc[ug.rc]));
	}
	if ((*ug.yfn)[ug.fId] >= 0)
	{
		iinv.Vdvv[ug.fId] = ((*ug.farea)[ug.fId] / (iinv.spc[ug.lc]));
	}
	else
	{
		iinv.Vdvv[ug.fId] = ((*ug.farea)[ug.fId] / (iinv.spc[ug.rc]));
	}
	if ((*ug.zfn)[ug.fId] >= 0)
	{
		iinv.Vdvw[ug.fId] = ((*ug.farea)[ug.fId] / (iinv.spc[ug.lc]));
	}
	else
	{
		iinv.Vdvw[ug.fId] = ((*ug.farea)[ug.fId] / (iinv.spc[ug.rc]));
	}
	
	iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] ;*/
	
	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.ajp[ug.fId] = iinv.rf[ug.fId] * (iinv.Vdvu[ug.fId] * (*ug.xfn)[ug.fId] * (*ug.xfn)[ug.fId] + iinv.Vdvv[ug.fId] * (*ug.yfn)[ug.fId] * (*ug.yfn)[ug.fId] + iinv.Vdvw[ug.fId] * (*ug.zfn)[ug.fId] * (*ug.zfn)[ug.fId]) * (*ug.farea)[ug.fId] / iinv.dist;
}

void INsInvterm::CalcINsBcFaceCorrectPresscoef()
{
	iinv.Vdvu[ug.fId] =  0;  // (VP / DV) J, used to calculate the surface velocity correction
	iinv.Vdvv[ug.fId] = 0;
	iinv.Vdvw[ug.fId] = 0;
	iinv.ajp[ug.fId] = 0;
}


  














EndNameSpace
