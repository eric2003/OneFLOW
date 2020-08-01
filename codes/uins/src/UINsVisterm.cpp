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

#include "UINsVisterm.h"
#include "INsInvterm.h"
#include "UINsInvterm.h"
#include "INsVisterm.h"
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


void UINsVisterm::CalcViscoff()
{
    if ( vis_model.vismodel == 0 ) return;
    ug.Init();
    uinsf.Init();
    visQ.Init( nscom.nEqu );
    //visT.Init( nscom.nTModel );
    Ivis.Init();
   // heat_flux.Init();

    Alloc();

    //this->SetVisPointer();

    this->PrepareField();
    this->CalcVisterm();
    //this->Addterm();

    DeAlloc();
}

void UINsVisterm::Alloc()
{
    visflux = new MRField( nscom.nEqu, ug.nFace );
}

void UINsVisterm::DeAlloc()
{
    delete visflux;
}

void UINsVisterm::PrepareField()
{
	uins_grad.Init();
	uins_grad.CalcGrad();  //计算梯度
    //ut_grad.CalcGradDebug();
}

void UINsVisterm::CalcVisterm()
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

        if ( ug.lc == 11 || ug.rc == 11 )
        {
        }

        //this->PrepareFaceValue();

        this->CalcFaceVisterm();  //要改动

    }
}

void UINsVisterm::CalcFaceVisterm()
{
    this->CalcNsVisterm();  //要改动
}


void UINsVisterm::CalcNsVisterm()
{

	Real l2rdx = (*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc];  //界面左右单元中心距
	Real l2rdy = (*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc];
	Real l2rdz = (*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc];


	Real Fn = (1/2) / (gcom.xfn * l2rdx + gcom.yfn * l2rdy + gcom.zfn * l2rdz) * gcom.farea;   // μ / ( n * d ) 法向扩散项系数
	Real Ft = (1/2) * ((visQ.dqdx[IIDX::IIU] * gcom.xfn + visQ.dqdy[IIDX::IIV] * gcom.yfn + visQ.dqdz[IIDX::IIW] * gcom.zfn) -
		(visQ.dqdx[IIDX::IIU] * l2rdx + visQ.dqdy[IIDX::IIV] * l2rdy + visQ.dqdz[IIDX::IIW] * l2rdz) / (gcom.xfn * l2rdx + gcom.yfn * l2rdy + gcom.zfn * l2rdz)) * gcom.farea;//归入源项的扩散项

	iinv.ai1[ug.lc] = iinv.ai1[ug.lc] + Fn;
	iinv.ai2[ug.rc] = iinv.ai2[ug.rc] + Fn;
	iinv.bm[ug.fId] = Ft;  //界面上归入源项的扩散项


	//iinv.flux[IIDX::IIRU] += (iinv.rm * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * half * (gcom.cvol1 + gcom.cvol2) / dist) * ((Pd1 + Pd2) - (iinv.pr - iinv.pl)) * gcom.farea;
	//iinv.flux[IIDX::IIRV] += (iinv.rm * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * half * (gcom.cvol1 + gcom.cvol2) / dist) * ((Pd1 + Pd2) - (iinv.pr - iinv.pl)) * gcom.farea;
	//iinv.flux[IIDX::IIRW] += (iinv.rm * SQR(gcom.xfn, gcom.yfn, gcom.zfn) * half * (gcom.cvol1 + gcom.cvol2) / dist) * ((Pd1 + Pd2) - (iinv.pr - iinv.pl)) * gcom.farea;

}

void UINsVisterm::CalcSrc()
{
	for (int cId = 0; cId < ug.nCell; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = ( *ug.lcf )[ug.fId];
			ug.rc = ( *ug.rcf )[ug.fId];

			iinv.buc[ug.cId] += iinv.bm[ug.fId];  //计入源项的扩散项
			iinv.bvc[ug.cId] += iinv.bm[ug.fId];  
			iinv.bwc[ug.cId] += iinv.bm[ug.fId];
			
			iinv.sp[ug.cId] = iinv.sp[ug.cId] + iinv.ai1[ug.lc];  //以cId单元构造动量方程时的系数（与质量通量相关）
			iinv.spj[ug.cId] = iinv.spj[ug.cId] + iinv.ai2[ug.rc]; //与cId相邻单元的系数
			iinv.sp2[ug.cId] = iinv.sp2[ug.rc] + iinv.ai2[ug.rc] + iinv.rr * gcom.cvol2 / (*uinsf.timestep)[0][ug.cId];
		}
		//iinv.spu[ug.cId] += iinv.rl * gcom.cvol / (*uinsf.timestep)[0][ug.cId]- gcom.cvol * visQ.dqdx1[IIDX::IIP]; //加上时间项和压力梯度项的cId单元系数(计算速度u)
		//iinv.spv[ug.cId] += iinv.rl * gcom.cvol / (*uinsf.timestep)[0][ug.cId] - gcom.cvol * visQ.dqdy1[IIDX::IIP];
		iinv.sp[ug.cId] += iinv.rl * gcom.cvol / (*uinsf.timestep)[0][ug.cId] - gcom.cvol * visQ.dqdz1[IIDX::IIP];
		iinv.buc[ug.cId] += iinv.rl * gcom.cvol * iinv.prim[IIDX::IIU] / (*uinsf.timestep)[0][ug.cId];
		iinv.bvc[ug.cId] += iinv.rl * gcom.cvol * iinv.prim[IIDX::IIV] / (*uinsf.timestep)[0][ug.cId];
		iinv.bwc[ug.cId] += iinv.rl * gcom.cvol * iinv.prim[IIDX::IIW] / (*uinsf.timestep)[0][ug.cId];
	}
}

//void UINsVisterm::AddVisFlux()
//{
 //   UnsGrid * grid = Zone::GetUnsGrid();
 //   MRField * res = GetFieldPointer< MRField >( grid, "res" );

 //   ONEFLOW::AddF2CField( res, visflux );
//}

void UINsVisterm::PrepareFaceValue()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];
    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.CalcTangent();

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.dqdx1[ iEqu ] = ( * uinsf.dqdx )[ iEqu ][ ug.lc ];
        visQ.dqdy1[ iEqu ] = ( * uinsf.dqdy )[ iEqu ][ ug.lc ];
        visQ.dqdz1[ iEqu ] = ( * uinsf.dqdz )[ iEqu ][ ug.lc ];

        visQ.dqdx2[ iEqu ] = ( * uinsf.dqdx )[ iEqu ][ ug.rc ];
        visQ.dqdy2[ iEqu ] = ( * uinsf.dqdy )[ iEqu ][ ug.rc ];
        visQ.dqdz2[ iEqu ] = ( * uinsf.dqdz )[ iEqu ][ ug.rc ];
    }

   // for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    //{
   //     visT.dqdx1[ iEqu ] = ( * uinsf.dtdx )[ iEqu ][ ug.lc ];
   //     visT.dqdy1[ iEqu ] = ( * uinsf.dtdy )[ iEqu ][ ug.lc ];
   //     visT.dqdz1[ iEqu ] = ( * uinsf.dtdz )[ iEqu ][ ug.lc ];

//        visT.dqdx2[ iEqu ] = ( * uinsf.dtdx )[ iEqu ][ ug.rc ];
//        visT.dqdy2[ iEqu ] = ( * uinsf.dtdy )[ iEqu ][ ug.rc ];
//        visT.dqdz2[ iEqu ] = ( * uinsf.dtdz )[ iEqu ][ ug.rc ];
//    }

    nscom.visl1 = ( * uinsf.visl )[ 0 ][ ug.lc ];
    nscom.visl2 = ( * uinsf.visl )[ 0 ][ ug.rc ];

    nscom.vist1 = ( * uinsf.vist )[ 0 ][ ug.lc ];
    nscom.vist2 = ( * uinsf.vist )[ 0 ][ ug.rc ];

    nscom.visl = half * ( nscom.visl1 + nscom.visl2 );
    nscom.vist = half * ( nscom.vist1 + nscom.vist2 );
    nscom.vis  = nscom.visl + nscom.vist;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.q1[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
        visQ.q2[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.rc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        visQ.q11[ iEqu ] = visQ.q1[ iEqu ];
        visQ.q22[ iEqu ] = visQ.q2[ iEqu ];
    }

   // for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
   // {
   //     visT.q1[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.lc ];
   //     visT.q2[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.rc ];
   // }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        visT.q11[ iEqu ] = visT.q1[ iEqu ];
        visT.q22[ iEqu ] = visT.q2[ iEqu ];
    }

    this->AverGrad();
    this->CalcFaceWeight();
    this->SaveFacePara();
}

void UINsVisterm::SaveFacePara()
{
    Ivis.dudx  = visQ.dqdx[ IIDX::IIU ];
    Ivis.dudy  = visQ.dqdy[ IIDX::IIU ];
    Ivis.dudz  = visQ.dqdz[ IIDX::IIU ];

    Ivis.dvdx  = visQ.dqdx[ IIDX::IIV ];
    Ivis.dvdy  = visQ.dqdy[ IIDX::IIV ];
    Ivis.dvdz  = visQ.dqdz[ IIDX::IIV ];

    Ivis.dwdx  = visQ.dqdx[ IIDX::IIW ];
    Ivis.dwdy  = visQ.dqdy[ IIDX::IIW ];
    Ivis.dwdz  = visQ.dqdz[ IIDX::IIW ];

	Ivis.dpdx = visQ.dqdx[IIDX::IIP];
	Ivis.dpdy = visQ.dqdy[IIDX::IIP];
	Ivis.dpdz = visQ.dqdz[IIDX::IIP];

	Ivis.p1 = visQ.q1[IIDX::IIP];
	Ivis.p2 = visQ.q2[IIDX::IIP];

    Ivis.um  = visQ.q[ IIDX::IIU ];
    Ivis.vm  = visQ.q[ IIDX::IIV ];
    Ivis.wm  = visQ.q[ IIDX::IIW ];

    //Ivis.dtdn = visT.dqdn[ IIDX::IITT ];
    //Ivis.tmid = visT.q[ IIDX::IITT ];
}

void UINsVisterm::CalcFaceWeight()
{
    vgg.CalcFaceWeight();
}


void UINsVisterm::CalcGradCoef()
{
    vgg.CalcGradCoef();
}


void UINsVisterm::PrepareCellGeom()
{
    vgg.PrepareCellGeom();
}


void ICalcLaminarViscosity(int flag)
{
	ug.Init();
	uinsf.Init();
	ug.SetStEd(flag);

	Real minLimit = 0.0;

	for (int cId = ug.ist; cId < ug.ied; ++cId)
	{
		//Real temperature = ( *uinsf.tempr )[ IIDX::IITT ][ cId ];
		//Real visl = Iutherland::ICalcViscosity( temperature );
		//( *uinsf.visl )[ 0 ][ cId ] = MAX( minLimit, visl );
	}
}

EndNameSpace

