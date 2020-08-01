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

#include "UINsVisterm.h"
#include "INsInvterm.h"
#include "UINsInvterm.h"
#include "INsVisterm.h"
#include "Iteration.h"
#include "HeatFlux.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "INsCtrl.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "VisGrad.h"
#include "UINsGrad.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "ULimiter.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


UINsVisterm::UINsVisterm()
{
	;
}

UINsVisterm::~UINsVisterm()
{
    ;
}


void UINsVisterm::CmpViscoff()
{
    if ( vis_model.vismodel == 0 ) return;
    ug.Init();
    uinsf.Init();
    visQ.Init( inscom.nEqu );

    //Alloc();

    this->PrepareField();
    this->CmpVisterm();

    //DeAlloc();
}

void UINsVisterm::Alloc()
{
    visflux = new MRField( inscom.nEqu, ug.nFace );
}

void UINsVisterm::DeAlloc()
{
    delete visflux;
}

void UINsVisterm::PrepareField()
{
	//uins_grad.Init();
	//uins_grad.CmpGrad();  //计算梯度
    //ut_grad.CmpGradDebug();
	this->CmpPreandVisGrad();
}

void UINsVisterm::CmpPreandVisGrad()
{
	(*uinsf.dqdx)[IIDX::IIR] = 0;
	(*uinsf.dqdy)[IIDX::IIR] = 0;
	(*uinsf.dqdz)[IIDX::IIR] = 0;

	(*uinsf.dqdx)[IIDX::IIU] = 0;
	(*uinsf.dqdy)[IIDX::IIU] = 0;
	(*uinsf.dqdz)[IIDX::IIU] = 0;

	(*uinsf.dqdx)[IIDX::IIV] = 0;
	(*uinsf.dqdy)[IIDX::IIV] = 0;
	(*uinsf.dqdz)[IIDX::IIV] = 0;

	(*uinsf.dqdx)[IIDX::IIW] = 0;
	(*uinsf.dqdy)[IIDX::IIW] = 0;
	(*uinsf.dqdz)[IIDX::IIW] = 0;

	(*uinsf.dqdx)[IIDX::IIP] = 0;
	(*uinsf.dqdy)[IIDX::IIP] = 0;
	(*uinsf.dqdz)[IIDX::IIP] = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		if (fId == 432)
		{
			int kkk = 1;
		}
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		Real dxl = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
		Real dyl = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
		Real dzl = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

		Real dxr = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
		Real dyr = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
		Real dzr = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

		Real delt1 = DIST(dxl, dyl, dzl);
		Real delt2 = DIST(dxr, dyr, dzr);
		Real delta = 1.0 / (delt1 + delt2);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;

		Real value1 = cl * (*uinsf.q)[IIDX::IIU][ug.lc] + cr * (*uinsf.q)[IIDX::IIU][ug.rc];
		Real value2 = cl * (*uinsf.q)[IIDX::IIV][ug.lc] + cr * (*uinsf.q)[IIDX::IIV][ug.rc];
		Real value3 = cl * (*uinsf.q)[IIDX::IIW][ug.lc] + cr * (*uinsf.q)[IIDX::IIW][ug.rc];
		Real value4 = cl * (*uinsf.q)[IIDX::IIP][ug.lc] + cr * (*uinsf.q)[IIDX::IIP][ug.rc];

		/*Real value1 = 0.5 * (*uinsf.q)[IIDX::IIU][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIU][ug.rc];
		Real value2 = 0.5 * (*uinsf.q)[IIDX::IIV][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIV][ug.rc];
		Real value3 = 0.5 * (*uinsf.q)[IIDX::IIW][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIW][ug.rc];
		Real value4 = 0.5 * (*uinsf.q)[IIDX::IIP][ug.lc] + 0.5 * (*uinsf.q)[IIDX::IIP][ug.rc];*/

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		(*uinsf.dqdx)[IIDX::IIU][ug.lc] += fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.lc] += fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.lc] += fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.lc] += fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.lc] += fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.lc] += fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.lc] += fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.lc] += fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.lc] += fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.lc] += fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.lc] += fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.lc] += fnza * value4;

		if (ug.fId < ug.nBFace) continue;
		(*uinsf.dqdx)[IIDX::IIU][ug.rc] += -fnxa * value1;
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] += -fnya * value1;
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] += -fnza * value1;
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] += -fnxa * value2;
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] += -fnya * value2;
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] += -fnza * value2;
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] += -fnxa * value3;
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] += -fnya * value3;
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] += -fnza * value3;
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] += -fnxa * value4;
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] += -fnya * value4;
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] += -fnza * value4;
	}

	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
		(*uinsf.dqdx)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIU][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIV][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIW][ug.cId] *= ovol;
		(*uinsf.dqdx)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdy)[IIDX::IIP][ug.cId] *= ovol;
		(*uinsf.dqdz)[IIDX::IIP][ug.cId] *= ovol;
	}

	for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];
		//if (ug.rc > ug.nCell)
		//{
		(*uinsf.dqdx)[IIDX::IIU][ug.rc] = (*uinsf.dqdx)[IIDX::IIU][ug.lc];
		(*uinsf.dqdy)[IIDX::IIU][ug.rc] = (*uinsf.dqdy)[IIDX::IIU][ug.lc];
		(*uinsf.dqdz)[IIDX::IIU][ug.rc] = (*uinsf.dqdz)[IIDX::IIU][ug.lc];
		(*uinsf.dqdx)[IIDX::IIV][ug.rc] = (*uinsf.dqdx)[IIDX::IIV][ug.lc];
		(*uinsf.dqdy)[IIDX::IIV][ug.rc] = (*uinsf.dqdy)[IIDX::IIV][ug.lc];
		(*uinsf.dqdz)[IIDX::IIV][ug.rc] = (*uinsf.dqdz)[IIDX::IIV][ug.lc];
		(*uinsf.dqdx)[IIDX::IIW][ug.rc] = (*uinsf.dqdx)[IIDX::IIW][ug.lc];
		(*uinsf.dqdy)[IIDX::IIW][ug.rc] = (*uinsf.dqdy)[IIDX::IIW][ug.lc];
		(*uinsf.dqdz)[IIDX::IIW][ug.rc] = (*uinsf.dqdz)[IIDX::IIW][ug.lc];
		(*uinsf.dqdx)[IIDX::IIP][ug.rc] = (*uinsf.dqdx)[IIDX::IIP][ug.lc];
		(*uinsf.dqdy)[IIDX::IIP][ug.rc] = (*uinsf.dqdy)[IIDX::IIP][ug.lc];
		(*uinsf.dqdz)[IIDX::IIP][ug.rc] = (*uinsf.dqdz)[IIDX::IIP][ug.lc];
		//}

	}

}


void UINsVisterm::CmpVisterm()
{
    for ( int fId = 0; fId < ug.nFace; ++ fId )
    {
        ug.fId = fId;

        ug.lc = ( * ug.lcf )[ ug.fId ];
        ug.rc = ( * ug.rcf )[ ug.fId ];

        if ( fId == 147489 )
        {
            int kkk = 1;
        }
		
		//iinv.ukl[ug.fId] = (*limf->qf1)[IIDX::IIU][ug.fId];
		//iinv.ukr[ug.fId] = (*limf->qf2)[IIDX::IIU][ug.fId];
		//iinv.vkl[ug.fId] = (*limf->qf1)[IIDX::IIV][ug.fId];
		//iinv.vkr[ug.fId] = (*limf->qf2)[IIDX::IIV][ug.fId];
		//iinv.wkl[ug.fId] = (*limf->qf1)[IIDX::IIW][ug.fId];
		//iinv.wkr[ug.fId] = (*limf->qf2)[IIDX::IIW][ug.fId];

        this->CmpFaceVisterm();  //要改动

    }

	/*for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		if (fId == 147489)
		{
			int kkk = 1;
		}

		//iinv.ukl[ug.fId] = (*limf->qf1)[IIDX::IIU][ug.fId];
		//iinv.ukr[ug.fId] = (*limf->qf2)[IIDX::IIU][ug.fId];
		//iinv.vkl[ug.fId] = (*limf->qf1)[IIDX::IIV][ug.fId];
		//iinv.vkr[ug.fId] = (*limf->qf2)[IIDX::IIV][ug.fId];
		//iinv.wkl[ug.fId] = (*limf->qf1)[IIDX::IIW][ug.fId];
		//iinv.wkr[ug.fId] = (*limf->qf2)[IIDX::IIW][ug.fId];

		this->CmpBcFaceVisterm();  //要改动

	}*/

}

void UINsVisterm::CmpFaceVisterm()
{

	iinv.l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	iinv.l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	iinv.l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];

	iinv.c2d = sqrt(iinv.l2rdx * iinv.l2rdx + iinv.l2rdy * iinv.l2rdy + iinv.l2rdz * iinv.l2rdz);

	iinv.vis = 1 / inscom.reynolds;  //动力粘度

	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.Fn[ug.fId] = iinv.vis * (*ug.farea)[ug.fId] / iinv.dist;

	Real dx1 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];
	Real dy1 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	Real dz1 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	Real dx2 = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.rc];
	Real dy2 = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.rc];
	Real dz2 = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.rc];

	Real de1 = DIST(dx1, dy1, dz1);
	Real de2 = DIST(dx2, dy2, dz2);
	Real de = 1.0 / (de1 + de2);

	iinv.f1[ug.fId] = de2 * de;  //左单元权重
	iinv.f2[ug.fId] = de1 * de;  //右单元权重

	iinv.Puf = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId];  //▽q*n

	iinv.Pvf = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pwf = (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*(*ug.xfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.yfn)[ug.fId] +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId];
	
	iinv.Pdu = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*iinv.l2rdx+
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*iinv.l2rdy+
	 	                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*iinv.l2rdz)/ iinv.dist;

	iinv.Pdv = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*iinv.l2rdx  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*iinv.l2rdy  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*iinv.l2rdz) / iinv.dist;

	iinv.Pdw = -((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.rc])*iinv.l2rdx  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*iinv.l2rdy  +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*iinv.l2rdz ) / iinv.dist;

	
	iinv.Ftu1 = iinv.Puf *(*ug.farea)[ug.fId]* iinv.vis;   //扩散项中归入源项的部分1
    iinv.Ftv1 = iinv.Pvf *(*ug.farea)[ug.fId]* iinv.vis;
	iinv.Ftw1 = iinv.Pwf *(*ug.farea)[ug.fId]* iinv.vis;

	iinv.Ftu2 = iinv.Pdu*(*ug.farea)[ug.fId] * iinv.vis;   //扩散项中归入源项的部分2
	iinv.Ftv2 = iinv.Pdv*(*ug.farea)[ug.fId] * iinv.vis;
	iinv.Ftw2 = iinv.Pdw*(*ug.farea)[ug.fId] * iinv.vis;



	iinv.PufT = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PvfT = ((iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);

	iinv.PwfT = ((iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIU][ug.rc])*(*ug.xfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIV][ug.rc])*(*ug.yfn)[ug.fId] +
		                 (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc])*(*ug.zfn)[ug.fId]);



	iinv.Pud = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.xfn)[ug.fId];

    iinv.Pvd = ((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		                (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.yfn)[ug.fId];

	iinv.Pwd =((iinv.f1[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdx)[IIDX::IIU][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdy)[IIDX::IIV][ug.rc]) +
		               (iinv.f1[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.lc] + iinv.f2[ug.fId] * (*uinsf.dqdz)[IIDX::IIW][ug.rc]))*(*ug.zfn)[ug.fId];

	iinv.FtuT = iinv.PufT * (*ug.farea)[ug.fId] * iinv.vis;  //Г(▽V)T，表面源项
	iinv.FtvT = iinv.PvfT * (*ug.farea)[ug.fId] * iinv.vis;
	iinv.FtwT = iinv.PwfT * (*ug.farea)[ug.fId] * iinv.vis;

	iinv.ai[ug.fId][0] += iinv.Fn[ug.fId];
	iinv.ai[ug.fId][1] += iinv.Fn[ug.fId];
	
	iinv.biu[ug.fId][0] = iinv.Ftu1 + iinv.Ftu2;
	iinv.biu[ug.fId][1] = -iinv.Ftu1 - iinv.Ftu2;

	iinv.biv[ug.fId][0] = iinv.Ftv1 + iinv.Ftv2;
	iinv.biv[ug.fId][1] = -iinv.Ftv1 - iinv.Ftv2;

	iinv.biw[ug.fId][0] = iinv.Ftw1 + iinv.Ftw2;
	iinv.biw[ug.fId][1] = -iinv.Ftw1 - iinv.Ftw2;
}

/*void UINsVisterm::CmpBcFaceVisterm()
{
	iinv.l2rdx[ug.fId] = (*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	iinv.l2rdy[ug.fId] = (*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc];
	iinv.l2rdz[ug.fId] = (*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc];

	iinv.c2d = sqrt(iinv.l2rdx[ug.fId] * iinv.l2rdx[ug.fId] + iinv.l2rdy[ug.fId] * iinv.l2rdy[ug.fId] + iinv.l2rdz[ug.fId] * iinv.l2rdz[ug.fId]);

	iinv.dist[ug.fId] = (*ug.xfn)[ug.fId] * ((*ug.xfc)[ug.fId] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.yfc)[ug.fId] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zfc)[ug.fId] - (*ug.zcc)[ug.lc]);

	iinv.Fn[ug.fId] = iinv.vis * (*ug.farea)[ug.fId] / 2*iinv.dist[ug.fId];


	iinv.Puf[ug.fId] = ((*uinsf.dqdx)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIU][ug.lc])*(*ug.yfn)[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIU][ug.lc])*(*ug.zfn)[ug.fId];  //▽q*n

	iinv.Pvf[ug.fId] = ((*uinsf.dqdx)[IIDX::IIV][ug.lc] )*(*ug.xfn)[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIV][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.Pwf[ug.fId] = ((*uinsf.dqdx)[IIDX::IIW][ug.lc])*(*ug.xfn)[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIW][ug.lc])*(*ug.yfn)[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId];

	iinv.Pdu[ug.fId] = -(((*uinsf.dqdx)[IIDX::IIU][ug.lc])*iinv.l2rdx[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIU][ug.lc] )*iinv.l2rdy[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIU][ug.lc] )*iinv.l2rdz[ug.fId]) / iinv.dist[ug.fId];

	iinv.Pdv[ug.fId] = -(((*uinsf.dqdx)[IIDX::IIV][ug.lc])*iinv.l2rdx[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIV][ug.lc])*iinv.l2rdy[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIV][ug.lc])*iinv.l2rdz[ug.fId]) / iinv.dist[ug.fId];

	iinv.Pdw[ug.fId] = -(((*uinsf.dqdx)[IIDX::IIW][ug.lc])*iinv.l2rdx[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIW][ug.lc])*iinv.l2rdy[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIW][ug.lc] )*iinv.l2rdz[ug.fId]) / iinv.dist[ug.fId];

	iinv.Ftu1[ug.fId] = iinv.Puf[ug.fId] * (*ug.farea)[ug.fId] * iinv.visu[ug.fId];   //扩散项中归入源项的部分1
	iinv.Ftv1[ug.fId] = iinv.Pvf[ug.fId] * (*ug.farea)[ug.fId] * iinv.visv[ug.fId];
	iinv.Ftw1[ug.fId] = iinv.Pwf[ug.fId] * (*ug.farea)[ug.fId] * iinv.visw[ug.fId];

	iinv.Ftu2[ug.fId] = iinv.Pdu[ug.fId] * (*ug.farea)[ug.fId] * iinv.visu[ug.fId];   //扩散项中归入源项的部分2
	iinv.Ftv2[ug.fId] = iinv.Pdv[ug.fId] * (*ug.farea)[ug.fId] * iinv.visv[ug.fId];
	iinv.Ftw2[ug.fId] = iinv.Pdw[ug.fId] * (*ug.farea)[ug.fId] * iinv.visw[ug.fId];



	iinv.PufT[ug.fId] = (((*uinsf.dqdx)[IIDX::IIU][ug.lc] )*(*ug.xfn)[ug.fId] +
		((*uinsf.dqdx)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		((*uinsf.dqdx)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId]);

	iinv.PvfT[ug.fId] = (((*uinsf.dqdy)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		((*uinsf.dqdy)[IIDX::IIW][ug.lc])*(*ug.zfn)[ug.fId]);

	iinv.PwfT[ug.fId] = (((*uinsf.dqdz)[IIDX::IIU][ug.lc])*(*ug.xfn)[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIV][ug.lc])*(*ug.yfn)[ug.fId] +
		((*uinsf.dqdz)[IIDX::IIW][ug.lc] )*(*ug.zfn)[ug.fId]);



	iinv.Pud[ug.fId] = (((*uinsf.dqdx)[IIDX::IIU][ug.lc]) +
		((*uinsf.dqdy)[IIDX::IIV][ug.lc]) +
		((*uinsf.dqdz)[IIDX::IIW][ug.lc]))*(*ug.xfn)[ug.fId];

	iinv.Pvd[ug.fId] = (((*uinsf.dqdx)[IIDX::IIU][ug.lc]) +
		((*uinsf.dqdy)[IIDX::IIV][ug.lc]) +
		((*uinsf.dqdz)[IIDX::IIW][ug.lc]))*(*ug.yfn)[ug.fId];

	iinv.Pwd[ug.fId] = (((*uinsf.dqdx)[IIDX::IIU][ug.lc]) +
		((*uinsf.dqdy)[IIDX::IIV][ug.lc]) +
		((*uinsf.dqdz)[IIDX::IIW][ug.lc]))*(*ug.zfn)[ug.fId];

	iinv.FtuT[ug.fId] = iinv.PufT[ug.fId] * (*ug.farea)[ug.fId] * iinv.visu[ug.fId];  //Г(▽V)T，表面源项
	iinv.FtvT[ug.fId] = iinv.PvfT[ug.fId] * (*ug.farea)[ug.fId] * iinv.visv[ug.fId];
	iinv.FtwT[ug.fId] = iinv.PwfT[ug.fId] * (*ug.farea)[ug.fId] * iinv.visw[ug.fId];

	iinv.ai[0][ug.fId] += iinv.Fn[ug.fId];
	iinv.ai[1][ug.fId] += iinv.Fn[ug.fId];

	iinv.biu[0][ug.fId] = iinv.Ftu1[ug.fId] + iinv.Ftu2[ug.fId];
	iinv.biu[1][ug.fId] = -iinv.Ftu1[ug.fId] - iinv.Ftu2[ug.fId];

	iinv.biv[0][ug.fId] = iinv.Ftv1[ug.fId] + iinv.Ftv2[ug.fId];
	iinv.biv[1][ug.fId] = -iinv.Ftv1[ug.fId] - iinv.Ftv2[ug.fId];

	iinv.biw[0][ug.fId] = iinv.Ftw1[ug.fId] + iinv.Ftw2[ug.fId];
	iinv.biw[1][ug.fId] = -iinv.Ftw1[ug.fId] - iinv.Ftw2[ug.fId];
}*/

void UINsVisterm::CmpUnsteadcoff()
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;

		iinv.spt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]/ iinv.timestep;  //矩阵对角线元素的非稳态项
		
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId]* iinv.wp[ug.cId]/ iinv.timestep;
		}
		else
		{
			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;
		}
	}

	/*for (int fId = 0; fId < ug.nBFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spt[ug.rc] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] / iinv.timestep;

		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.up[ug.rc] = ((*uinsf.q)[IIDX::IIU][ug.rc]+ (*uinsf.q)[IIDX::IIU][ug.rc])/2;
			iinv.vp[ug.rc] = ((*uinsf.q)[IIDX::IIV][ug.rc]+ (*uinsf.q)[IIDX::IIU][ug.rc])/2;
			iinv.wp[ug.rc] = ((*uinsf.q)[IIDX::IIW][ug.rc]+ (*uinsf.q)[IIDX::IIU][ug.rc])/2;


			iinv.but[ug.rc] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] * iinv.up[ug.rc] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.rc] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] * iinv.vp[ug.rc] / iinv.timestep;
			iinv.bwt[ug.rc] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] * iinv.wp[ug.rc] / iinv.timestep;

		}
		else
		{
			iinv.but[ug.cId] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] * iinv.up[ug.rc] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] * iinv.vp[ug.rc] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.farea)[ug.fId] * (*uinsf.q)[IIDX::IIR][ug.rc] * iinv.wp[ug.rc] / iinv.timestep;
		}

	}*/


	for (int cId = ug.nCell; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.spt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] / iinv.timestep;  //矩阵对角线元素的非稳态项
		
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
	
		{
			iinv.up[ug.cId] = (*uinsf.q)[IIDX::IIU][ug.cId];
			iinv.vp[ug.cId] = (*uinsf.q)[IIDX::IIV][ug.cId];
			iinv.wp[ug.cId] = (*uinsf.q)[IIDX::IIW][ug.cId];

			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;


		}
		else
		{
			iinv.but[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.up[ug.cId] / iinv.timestep; //源项的非稳态项
			iinv.bvt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.vp[ug.cId] / iinv.timestep;
			iinv.bwt[ug.cId] = (*ug.cvol)[ug.cId] * (*uinsf.q)[IIDX::IIR][ug.cId] * iinv.wp[ug.cId] / iinv.timestep;
		}

	}

}



void UINsVisterm::CmpINsSrc()
{
	iinv.spc = 0;
	iinv.buc = 0;
	iinv.bvc = 0;
	iinv.bwc = 0;

	for (int fId = 0; fId < ug.nFace; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spc[ug.lc] += iinv.ai[ug.fId][0];
		iinv.spc[ug.rc] += iinv.ai[ug.fId][1];

		iinv.buc[ug.lc] += iinv.biu[ug.fId][0];
		iinv.buc[ug.rc] += iinv.biu[ug.fId][1];

		iinv.bvc[ug.lc] += iinv.biv[ug.fId][0];
		iinv.bvc[ug.rc] += iinv.biv[ug.fId][1];

		iinv.bwc[ug.lc] += iinv.biw[ug.fId][0];
		iinv.bwc[ug.rc] += iinv.biw[ug.fId][1];

	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.spc[ug.cId] += iinv.spt[ug.cId];

		iinv.buc[ug.cId] += iinv.but[ug.cId]- (*ug.cvol)[ug.cId] * (*uinsf.dqdx)[IIDX::IIP][ug.cId];
		iinv.bvc[ug.cId] += iinv.bvt[ug.cId] -(*ug.cvol)[ug.cId] * (*uinsf.dqdy)[IIDX::IIP][ug.cId];
		iinv.bwc[ug.cId] += iinv.bwt[ug.cId] - (*ug.cvol)[ug.cId] * (*uinsf.dqdz)[IIDX::IIP][ug.cId];

		//cout << "iinv.buc=" << iinv.buc[ug.cId] <<"cId="<< ug.cId<< "\n";

		int fn = (*ug.c2f)[ug.cId].size();
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.sj.resize(ug.nTCell, fn);
			iinv.sd.resize(ug.nTCell, fn);
		}
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (ug.cId == ug.lc)
			{
				iinv.sj[ug.cId][iFace] = -iinv.ai[ug.fId][0];  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
				iinv.sd[ug.cId][iFace] = ug.rc;
			}
			else if (ug.cId == ug.rc)
			{
				iinv.sj[ug.cId][iFace] = -iinv.ai[ug.fId][1];  //矩阵非零系数，动量方程中与主单元相邻的单元面通量
				iinv.sd[ug.cId][iFace] = ug.lc;
			}

		}
	}
}



EndNameSpace

