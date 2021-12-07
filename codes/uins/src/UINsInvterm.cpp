/*---------------------------------------------------------------------------*\
	OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
	Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "UINsInvterm.h"
#include "INsInvterm.h"
#include "UINsVisterm.h"
#include "UINsGrad.h"
#include "BcData.h"
#include "INsBcSolver.h"
#include "Zone.h"
#include "Atmosphere.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "UCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "HXMath.h"
#include "Multigrid.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "UINsLimiter.h"
#include "FieldImp.h"
#include "Iteration.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "Ctrl.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace(ONEFLOW)

UINsInvterm::UINsInvterm()
{
	limiter = new INsLimiter();
	limf = limiter->limf;
}

UINsInvterm::~UINsInvterm()
{
	delete limiter;
}

void UINsInvterm::CalcLimiter()
{
	limiter->CalcLimiter();
}

void UINsInvterm::CalcInvFace()  //Cell data reconstruction
{
	//uins_grad.Init();
	//uins_grad.CalcGrad();


	this->CalcLimiter();   //Don't change it

	this->GetQlQrField();  //Don't change it

	//this->ReconstructFaceValueField();  //Don't change it

	this->BoundaryQlQrFixField();  //Don't change it
}

void UINsInvterm::GetQlQrField()
{
	limf->GetQlQr();
}

void UINsInvterm::ReconstructFaceValueField()
{
	limf->CalcFaceValue();
	//limf->CalcFaceValueWeighted();
}

void UINsInvterm::BoundaryQlQrFixField()
{
	limf->BcQlQrFix();
}

void UINsInvterm::CalcInvcoff()
{
	if (nscom.icmpInv == 0) return;
	iinv.Init();
	ug.Init();
	uinsf.Init();
	//Alloc();

	//this->CalcInvFace();
	this->CalcInvMassFlux();  //needs to be changed

   //DeAlloc();
}

void UINsInvterm::CalcINsTimeStep()
{
	iinv.timestep = GetDataValue< Real >("global_dt");
}

void UINsInvterm::CalcINsPreflux()
{
	if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
	//if (ctrl.currTime == 0.001 && Iteration::outerSteps == 1)
	{
		if (nscom.icmpInv == 0) return;
		iinv.Init();
		ug.Init();
		uinsf.Init();
		//Alloc();

		this->CalcInvFace();
		this->INsPreflux();

		//DeAlloc();
	}

}

void UINsInvterm::INsPreflux()
{
	this->Initflux();

	/*RealField *rf = new RealField(ug.nFaces);
	RealField *uf = new RealField(ug.nFaces);
	RealField *vf = new RealField(ug.nFaces);
	RealField *wf = new RealField(ug.nFaces);
	RealField *fq = new RealField(ug.nFaces);*/

	for (int fId = ug.nBFaces; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CalcINsinvFlux();

	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareFaceValue();

		this->CalcINsBcinvFlux();
	}

	/*delete iinv.rf;
	delete iinv.uf;
	delete iinv.vf;
	delete iinv.wf;
	delete iinv.fq;*/

}
void UINsInvterm::Initflux()
{
	iinv.f1.resize(ug.nFaces);
	iinv.f2.resize(ug.nFaces);
	iinv.rf.resize(ug.nFaces);
	iinv.uf.resize(ug.nFaces);
	iinv.vf.resize(ug.nFaces);
	iinv.wf.resize(ug.nFaces);
	iinv.Vdvu.resize(ug.nFaces);
	iinv.Vdvv.resize(ug.nFaces);
	iinv.Vdvw.resize(ug.nFaces);
	iinv.aju.resize(ug.nFaces);
	iinv.ajv.resize(ug.nFaces);
	iinv.ajw.resize(ug.nFaces);
	iinv.VdU.resize(ug.nTCell);
	iinv.VdV.resize(ug.nTCell);
	iinv.VdW.resize(ug.nTCell);
	iinv.buc.resize(ug.nTCell);
	iinv.bvc.resize(ug.nTCell);
	iinv.bwc.resize(ug.nTCell);
	iinv.bp.resize(ug.nTCell);
	iinv.ajp.resize(ug.nFaces);
	iinv.sju.resize(ug.nTCell);
	iinv.sjv.resize(ug.nTCell);
	iinv.sjw.resize(ug.nTCell);
	iinv.fq.resize(ug.nFaces);
	iinv.spc.resize(ug.nTCell);
	iinv.ai.resize(ug.nFaces,2);
	iinv.biu.resize(ug.nFaces,2);
	iinv.biv.resize(ug.nFaces,2);
	iinv.biw.resize(ug.nFaces,2);
	//iinv.sj.resize(ug.nTCell, 4);
	//iinv.sd.resize(ug.nTCell, 4);
	//iinv.sjp.resize(ug.nTCell, 4);
	//iinv.sjd.resize(ug.nTCell, 4);
	iinv.spp.resize(ug.nTCell);
	iinv.pp.resize(ug.nTCell);
	iinv.uu.resize(ug.nTCell);
	iinv.vv.resize(ug.nTCell);
	iinv.ww.resize(ug.nTCell);
	iinv.uuj.resize(ug.nFaces);
	iinv.vvj.resize(ug.nFaces);
	iinv.wwj.resize(ug.nFaces);
	iinv.muc.resize(ug.nTCell);
	iinv.mvc.resize(ug.nTCell);
	iinv.mwc.resize(ug.nTCell);
	iinv.mp.resize(ug.nTCell);
	iinv.uc.resize(ug.nTCell);
	iinv.vc.resize(ug.nTCell);
	iinv.wc.resize(ug.nTCell);
	iinv.up.resize(ug.nTCell);
	iinv.vp.resize(ug.nTCell);
	iinv.wp.resize(ug.nTCell);
	iinv.spt.resize(ug.nTCell);
	iinv.but.resize(ug.nTCell);
	iinv.bvt.resize(ug.nTCell);
	iinv.bwt.resize(ug.nTCell);
	iinv.dqqdx.resize(ug.nTCell);
	iinv.dqqdy.resize(ug.nTCell);
	iinv.dqqdz.resize(ug.nTCell);
	iinv.Fn.resize(ug.nFaces);
	iinv.Fnu.resize(ug.nFaces);
	iinv.Fnv.resize(ug.nFaces);
	iinv.Fnw.resize(ug.nFaces);
	iinv.Fpu.resize(ug.nFaces);
	iinv.Fpv.resize(ug.nFaces);
	iinv.Fpw.resize(ug.nFaces);
	iinv.dsrl.resize(ug.nFaces);
	iinv.elrn.resize(ug.nFaces);
	//iinv.value.resize(ug.nFaces);
	iinv.mu.resize(ug.nCells);
	iinv.mv.resize(ug.nCells);
	iinv.mw.resize(ug.nCells);
	iinv.mua.resize(ug.nCells);
	iinv.mva.resize(ug.nCells);
	iinv.mwa.resize(ug.nCells);
	iinv.res_pp.resize(ug.nCells);
	iinv.res_up.resize(ug.nCells);
	iinv.res_vp.resize(ug.nCells);
	iinv.res_wp.resize(ug.nCells);

	iinv.ai1 = 0;
	iinv.ai2 = 0;
	iinv.spu1 = 1;
	iinv.spv1 = 1;
	iinv.spw1 = 1;
	iinv.spu2 = 1;
	iinv.spv2 = 1;
	iinv.spw2 = 1;

	iinv.buc = 0;
	iinv.bvc = 0;
	iinv.bwc = 0;
	iinv.sp1 = 0;
	iinv.sp2 = 0;
	iinv.spj = 0;
	iinv.spp = 0;
	iinv.sppu = 0;
	iinv.sppv = 0;
	iinv.sppw = 0;

	iinv.bpu = 0;
	iinv.bpv = 0;
	iinv.bpw = 0;
	iinv.bp = 0;
	iinv.pp = 0;

	iinv.muc = 0;
	iinv.mvc = 0;
	iinv.mwc = 0;

	iinv.uc = 0;
	iinv.vc = 0;
	iinv.wc = 0;

	iinv.uu = 0;
	iinv.vv = 0;
	iinv.ww = 0;
}

void UINsInvterm::CalcInvMassFlux()
{

	for (int fId = 0; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CalcINsinvTerm();
	}
}

void UINsInvterm::PrepareFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];

	nscom.gama1 = ( * uinsf.gama )[ 0 ][ ug.lc ];
	nscom.gama2 = ( * uinsf.gama )[ 0 ][ ug.rc ];
	nscom.gama  = half * ( nscom.gama1 + nscom.gama2 );

	iinv.gama1 = nscom.gama1;
	iinv.gama2 = nscom.gama2;

	for (int iEqu = 0; iEqu < limf->nEqu; ++iEqu)
	{
		iinv.prim1[iEqu] = (*limf->q)[iEqu][ug.lc];
		iinv.prim2[iEqu] = (*limf->q)[iEqu][ug.rc];
	}
}

void UINsInvterm::PrepareProFaceValue()
{
	gcom.xfn = (*ug.xfn)[ug.fId];
	gcom.yfn = (*ug.yfn)[ug.fId];
	gcom.zfn = (*ug.zfn)[ug.fId];
	gcom.vfn = (*ug.vfn)[ug.fId];
	gcom.farea = (*ug.farea)[ug.fId];


	//for (int iEqu = 0; iEqu < limf->nEqu; ++iEqu)
	//{
	iinv.prim1[IIDX::IIR] = (*uinsf.q)[IIDX::IIR][ug.lc];
	iinv.prim1[IIDX::IIU] = iinv.uc[ug.lc];
	iinv.prim1[IIDX::IIV] = iinv.vc[ug.lc];
	iinv.prim1[IIDX::IIW] = iinv.wc[ug.lc];
	iinv.prim1[IIDX::IIP] = (*uinsf.q)[IIDX::IIP][ug.lc];

	iinv.prim2[IIDX::IIR] = (*uinsf.q)[IIDX::IIR][ug.rc];
	iinv.prim2[IIDX::IIU] = iinv.uc[ug.rc];
	iinv.prim2[IIDX::IIV] = iinv.vc[ug.rc];
	iinv.prim2[IIDX::IIW] = iinv.vc[ug.rc];
	iinv.prim2[IIDX::IIP] = (*uinsf.q)[IIDX::IIP][ug.rc];
	//}
}

UINsInvterm NonZero;
void UINsInvterm::Init()
{
	this->Number = 0;
}

void UINsInvterm::MomPre()
{
	this->CalcINsMomRes();

	/*iinv.muc = 0;
	iinv.mvc = 0;
	iinv.mwc = 0;

	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;
		int fn = (*ug.c2f)[ug.cId].size();
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (ug.cId == ug.lc)
			{
				iinv.muc[ug.cId] += -iinv.sj[ug.cId][iFace] * (*uinsf.q)[IIDX::IIU][ug.rc];   //When Gauss Seidel iteration is used, the influence of adjacent elements on it is unnecessary for matrix method
				iinv.mvc[ug.cId] += -iinv.sj[ug.cId][iFace] * (*uinsf.q)[IIDX::IIV][ug.rc];
				iinv.mwc[ug.cId] += -iinv.sj[ug.cId][iFace] * (*uinsf.q)[IIDX::IIW][ug.rc];

			}
			else if (ug.cId == ug.rc)
			{

				iinv.muc[ug.cId] += -iinv.sj[ug.cId][iFace] * (*uinsf.q)[IIDX::IIU][ug.lc]; //When Gauss Seidel iteration is used, the influence of adjacent elements on it is unnecessary for matrix method
				iinv.mvc[ug.cId] += -iinv.sj[ug.cId][iFace] * (*uinsf.q)[IIDX::IIV][ug.lc];
				iinv.mwc[ug.cId] += -iinv.sj[ug.cId][iFace] * (*uinsf.q)[IIDX::IIW][ug.lc];

			}
		}


			iinv.uc[ug.cId] = (iinv.muc[ug.cId] + iinv.buc[ug.cId]) / (iinv.spc[ug.cId]);  //Predicted value of speed at the next moment


			iinv.vc[ug.cId] = (iinv.mvc[ug.cId] + iinv.bvc[ug.cId]) / (iinv.spc[ug.cId]);

			iinv.wc[ug.cId] = (iinv.mwc[ug.cId] + iinv.bwc[ug.cId]) / (iinv.spc[ug.cId]);

	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		BcInfo * bcInfo = ug.bcRecord->bcInfo;

		ug.fId = bcInfo->bcFace[ug.ir][fId];
		ug.bcNameId = bcInfo->bcNameId[ug.ir][fId];

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		nscom.bcdtkey = 0;
		if (ug.bcNameId == -1) return; //interface
		int dd = ns_bc_data.r2d[ug.bcNameId];
		if (dd != -1)
		{
			nscom.bcdtkey = 1;
			nscom.bcflow = &ns_bc_data.dataList[dd];
		}

		if (nscom.bcdtkey == 0)
		{
			iinv.uc[ug.rc] = -iinv.uc[ug.lc] + 2 * gcom.vfx;
			iinv.vc[ug.rc] = -iinv.vc[ug.lc] + 2 * gcom.vfy;
			iinv.wc[ug.rc] = -iinv.wc[ug.lc] + 2 * gcom.vfz;
		}
		else
		{
			iinv.uc[ug.rc] = -iinv.uc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIU];
			iinv.vc[ug.rc] = -iinv.vc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIV];
			iinv.wc[ug.rc] = -iinv.wc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIW];
		}
	}*/

	//Bgmres solution
	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{                                          
		int fn = (*ug.c2f)[cId].size();                             // Number of adjacent cells                             
		NonZero.Number += fn;                                       // The number of nonzero elements on a non diagonal line
	}
	NonZero.Number = NonZero.Number + ug.nTCell;                     // The total number of nonzero elements     
	Rank.RANKNUMBER = ug.nTCell;                                     // Row and column size of matrix
	Rank.NUMBER = NonZero.Number;                                    // The number of non-zero elements of matrix is transferred to the calculation program
	Rank.COLNUMBER = 1;                                              // Number of right end items
	Rank.Init();                                                     //Intermediate variable passed into GMRES calculation program
	double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + fn + 1;                                                  // The number of non-zero elements in the first n + 1 row
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // Number of adjacent faces
			ug.lc = (*ug.lcf)[fId];                                                                     // Face left unit
			ug.rc = (*ug.rcf)[fId];                                                                     // Face right unit
			if (cId == ug.lc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[fId][0];
				Rank.TempJA[n + iFace] = ug.rc;
			}
			else if (cId == ug.rc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[fId][1];
				Rank.TempJA[n + iFace] = ug.lc;
			}
		}
		Rank.TempA[n + fn] = iinv.spc[cId];                          //Primary diagonal element value
		Rank.TempJA[n + fn] = cId;                                      //Principal diagonal ordinate

	}
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.buc[cId];
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		iinv.uc[cId] = Rank.TempX[cId][0];                       // Output of solution
	}
	residual_u = Rank.residual;
	iinv.res_u = residual_u;

	Rank.Deallocate();
	//cout << "residual_u:" << residual_u << std::endl;




	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		int fn = (*ug.c2f)[cId].size();                             //Number of adjacent cells                                    
		NonZero.Number += fn;                                          //The number of nonzero elements on a non diagonal line
	}
	NonZero.Number = NonZero.Number + ug.nTCell;                     //The total number of nonzero elements         
	Rank.RANKNUMBER = ug.nTCell;                                     // Row and column size of matrix
	Rank.NUMBER = NonZero.Number;                                    // The number of non-zero elements of matrix is transferred to the calculation program
	Rank.COLNUMBER = 1;                                              //Number of right end items
	Rank.Init();                                                     //Intermediate variable passed into GMRES calculation program
	//double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + fn + 1;                                                  // The number of non-zero elements in the first n + 1 row
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // Number of adjacent faces
			ug.lc = (*ug.lcf)[fId];                                                                     // Face left unit
			ug.rc = (*ug.rcf)[fId];                                                                     // Face right unit
			if (cId == ug.lc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[fId][0];
				Rank.TempJA[n + iFace] = ug.rc;
			}
			else if (cId == ug.rc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[fId][1];
				Rank.TempJA[n + iFace] = ug.lc;
			}
		}
		Rank.TempA[n + fn] = iinv.spc[cId];                          //Primary diagonal element value
		Rank.TempJA[n + fn] = cId;                                      //Principal diagonal ordinate

	}


	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.bvc[cId];
	}	
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		iinv.vc[cId] = Rank.TempX[cId][0];
	}
	residual_v = Rank.residual;
	iinv.res_v = residual_v;

	Rank.Deallocate();

	//cout << "residual_v:" << residual_v << std::endl;


	NonZero.Number = 0;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		int fn = (*ug.c2f)[cId].size();                             //Number of adjacent cells                                    
		NonZero.Number += fn;                                          //The number of nonzero elements on a non diagonal line
	}
	NonZero.Number = NonZero.Number + ug.nTCell;                     //The total number of nonzero elements         
	Rank.RANKNUMBER = ug.nTCell;                                     // Row and column size of matrix
	Rank.NUMBER = NonZero.Number;                                    // The number of non-zero elements of matrix is transferred to the calculation program
	Rank.COLNUMBER = 1;                                              //Number of right end items
	Rank.Init();                                                     //Intermediate variable passed into GMRES calculation program
	//double residual_u, residual_v, residual_w;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + fn + 1;                                                  // The number of non-zero elements in the first n + 1 row
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                                                            // Number of adjacent faces
			ug.lc = (*ug.lcf)[fId];                                                                     // Face left unit
			ug.rc = (*ug.rcf)[fId];                                                                     // Face right unit
			if (cId == ug.lc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[fId][0];
				Rank.TempJA[n + iFace] = ug.rc;
			}
			else if (cId == ug.rc)
			{
				Rank.TempA[n + iFace] = -iinv.ai[fId][1];
				Rank.TempJA[n + iFace] = ug.lc;
			}
		}
		Rank.TempA[n + fn] = iinv.spc[cId];                          //Primary diagonal element value
		Rank.TempJA[n + fn] = cId;                                      //Principal diagonal ordinate

	}

	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		Rank.TempB[cId][0] = iinv.bwc[cId];
	}
	bgx.BGMRES();
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		iinv.wc[cId] = Rank.TempX[cId][0];
	}
	residual_w = Rank.residual;
	iinv.res_w = residual_w;

	Rank.Deallocate();

	//cout << "residual_w:" << residual_w << std::endl;

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		BcInfo * bcInfo = ug.bcRecord->bcInfo;

		ug.fId = bcInfo->bcFace[ ug.ir ][ fId ];
		ug.bcNameId = bcInfo->bcNameId[ ug.ir ][ fId ];

		ug.lc = ( * ug.lcf )[ ug.fId ];
		ug.rc = ( * ug.rcf )[ ug.fId ];

		if (ug.bctype < 0)
		{
			false;
		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
			nscom.bcdtkey = 0;
			if (ug.bcNameId == -1) return; //interface
			int dd = ins_bc_data.r2d[ug.bcNameId];
			if (dd != -1)
			{
				nscom.bcdtkey = 1;
				nscom.bcflow = &ins_bc_data.dataList[dd];
			}

			if (nscom.bcdtkey == 0)
			{
				iinv.uc[ug.rc] = -iinv.uc[ug.lc] + 2 * gcom.vfx;
				iinv.vc[ug.rc] = -iinv.vc[ug.lc] + 2 * gcom.vfy;
				iinv.wc[ug.rc] = -iinv.wc[ug.lc] + 2 * gcom.vfz;
			}
			else
			{
				iinv.uc[ug.rc] = -iinv.uc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIU];
				iinv.vc[ug.rc] = -iinv.vc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIV];
				iinv.wc[ug.rc] = -iinv.wc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIW];
			}

		}

		else if (ug.bctype == BC::INFLOW)
		{
			iinv.uc[ug.rc] = nscom.inflow[IIDX::IIU];
			iinv.vc[ug.rc] = nscom.inflow[IIDX::IIV];
			iinv.wc[ug.rc] = nscom.inflow[IIDX::IIW];
		}

		else if (ug.bctype == BC::OUTFLOW)
		{
			iinv.uc[ug.rc] = iinv.uc[ug.lc];
			iinv.vc[ug.rc] = iinv.vc[ug.lc];
			iinv.wc[ug.rc] = iinv.wc[ug.lc];
		}

		else if (ug.bctype == BC::POLE || ug.bctype / 10 == BC::POLE)
		{
			iinv.uc[ug.rc] = iinv.uc[ug.lc];
			iinv.vc[ug.rc] = iinv.vc[ug.lc];
			iinv.wc[ug.rc] = iinv.wc[ug.lc];
		}

		else if (ug.bctype == BC::EXTRAPOLATION)
		{
			iinv.uc[ug.rc] = iinv.uc[ug.lc];
			iinv.vc[ug.rc] = iinv.vc[ug.lc];
			iinv.wc[ug.rc] = iinv.wc[ug.lc];
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			Real vx1 = iinv.uc[ug.lc];
			Real vy1 = iinv.vc[ug.lc];
			Real vz1 = iinv.wc[ug.lc];

			Real vnRelative1 = (*ug.xfn)[ug.fId] * vx1 + (*ug.yfn)[ug.fId] * vy1 + (*ug.zfn)[ug.fId] * vz1 - (*ug.vfn)[ug.fId];

			iinv.uc[ug.rc] = iinv.uc[ug.lc]-two* (*ug.xfn)[ug.fId] * vnRelative1;
			iinv.vc[ug.rc] = iinv.vc[ug.lc]- two * (*ug.yfn)[ug.fId] * vnRelative1;
			iinv.wc[ug.rc] = iinv.wc[ug.lc]- two * (*ug.zfn)[ug.fId] * vnRelative1;

		}
		
		else if (ug.bctype == BC::FARFIELD)
		{
			Real rin = (*uinsf.q)[IIDX::IIR][ug.lc];
			Real uin = iinv.uc[ug.lc];
			Real vin = iinv.vc[ug.lc];
			Real win = iinv.wc[ug.lc];
			Real pin = (*uinsf.q)[IIDX::IIP][ug.lc];

			gcom.xfn *= nscom.faceOuterNormal;
			gcom.yfn *= nscom.faceOuterNormal;
			gcom.zfn *= nscom.faceOuterNormal;

			Real rref = nscom.inflow[IIDX::IIR];
			Real uref = nscom.inflow[IIDX::IIU];
			Real vref = nscom.inflow[IIDX::IIV];
			Real wref = nscom.inflow[IIDX::IIW];
			Real pref = nscom.inflow[IIDX::IIP];

			Real vnref = gcom.xfn * uref + gcom.yfn * vref + gcom.zfn * wref - gcom.vfn;
			Real vnin = gcom.xfn * uin + gcom.yfn * vin + gcom.zfn * win - (*ug.vfn)[ug.fId];

			Real cref = sqrt(ABS(nscom.gama_ref * pref / rref));
			Real cin = sqrt(ABS(nscom.gama * pin / rin));

			Real gamm1 = nscom.gama - one;

			Real velin = DIST(uin, vin, win);

			//Supersonic
			if (velin > cin)
			{
				if (vnin >= 0.0)
				{
					iinv.uc[ug.rc] = iinv.uc[ug.lc];
					iinv.vc[ug.rc] = iinv.vc[ug.lc];
					iinv.wc[ug.rc] = iinv.wc[ug.lc];
				}
				else
				{
					iinv.uc[ug.rc] = nscom.inflow[IIDX::IIU];
					iinv.vc[ug.rc] = nscom.inflow[IIDX::IIV];
					iinv.wc[ug.rc] = nscom.inflow[IIDX::IIW];
				}
			}
			else
			{
				//subsonic
				Real riemp = vnin + 2.0 * cin / gamm1;
				Real riemm = vnref - 2.0 * cref / gamm1;
				Real vnb = half * (riemp + riemm);
				Real cb = fourth * (riemp - riemm) * gamm1;

				Real vtx, vty, vtz, entr;
				if (vnb >= 0.0)
				{
					// exit
					entr = pin / pow(rin, nscom.gama);

					vtx = uin - gcom.xfn * vnin;
					vty = vin - gcom.yfn * vnin;
					vtz = win - gcom.zfn * vnin;
				}
				else
				{
					//inlet
					entr = pref / pow(rref, nscom.gama);
					vtx = uref - gcom.xfn * vnref;
					vty = vref - gcom.yfn * vnref;
					vtz = wref - gcom.zfn * vnref;
				}

				Real rb = pow((cb * cb / (entr * nscom.gama)), one / gamm1);
				Real ub = vtx + gcom.xfn * vnb;
				Real vb = vty + gcom.yfn * vnb;
				Real wb = vtz + gcom.zfn * vnb;
				Real pb = cb * cb * rb / nscom.gama;

				
				iinv.uc[ug.rc] = ub;
				iinv.vc[ug.rc] = vb;
				iinv.wc[ug.rc] = wb;

			}
		}

		else if (ug.bctype == BC::OVERSET)
		{
			;
		}

		else if (ug.bctype  == BC::GENERIC_2)
		{
		
			;
			
		}
		

	}

/*for (int cId = 0; cId < ug.nCells; cId++)
{
	ug.cId = cId;

	iinv.uc[ug.cId] = 0.0001;
	iinv.vc[ug.cId] = 0;
	iinv.wc[ug.cId] = 0;

}

for (int fId = 0; fId < ug.nBFaces; ++fId)
{
	ug.fId = fId;

	BcInfo * bcInfo = ug.bcRecord->bcInfo;

	ug.fId = bcInfo->bcFace[ug.ir][fId];
	ug.bcNameId = bcInfo->bcNameId[ug.ir][fId];

	ug.lc = (*ug.lcf)[ug.fId];
	ug.rc = (*ug.rcf)[ug.fId];

	nscom.bcdtkey = 0;
	if (ug.bcNameId == -1) return; //interface
	int dd = ns_bc_data.r2d[ug.bcNameId];
	if (dd != -1)
	{
		nscom.bcdtkey = 1;
		nscom.bcflow = &ns_bc_data.dataList[dd];
	}

	if (nscom.bcdtkey == 0)
	{
		iinv.uc[ug.rc] = -iinv.uc[ug.lc] + 2 * gcom.vfx;
		iinv.vc[ug.rc] = -iinv.vc[ug.lc] + 2 * gcom.vfy;
		iinv.wc[ug.rc] = -iinv.wc[ug.lc] + 2 * gcom.vfz;
	}
	else
	{
		iinv.uc[ug.rc] = -iinv.uc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIU];
		iinv.vc[ug.rc] = -iinv.vc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIV];
		iinv.wc[ug.rc] = -iinv.wc[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIW];
	}
}*/

	/*Output the residuals to a TXT file*/
	/*ofstream fileres_u("residual_u.txt", std::ios::app);
	//fileres_u << "residual_u:" << residual_u << std::endl;
	fileres_u << residual_u << std::endl;
	fileres_u.close();

	ofstream fileres_v("residual_v.txt", std::ios::app);
	//fileres_v << "residual_v:" << residual_v << std::endl;
	fileres_v << residual_v << std::endl;
	fileres_v.close();

	ofstream fileres_w("residual_w.txt", std::ios::app);
	//fileres_w << "residual_w:" << residual_w << std::endl;
	fileres_w <<residual_w << std::endl;
	fileres_w.close();*/
}

void UINsInvterm::CalcFaceflux()
{

	iinv.Init();
	ug.Init();
	uinsf.Init();
	//Alloc();
	//this->CalcInvFace();  //Boundary treatment
	for (int fId = ug.nBFaces; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareProFaceValue();

		this->CalcINsFaceflux();
	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->PrepareProFaceValue();

		this->CalcINsBcFaceflux();
	}

}

void UINsInvterm::CalcINsMomRes()
{
	iinv.res_u = 0;
	iinv.res_v = 0;
	iinv.res_w = 0;

	//Conditions for judging convergence of iteration
	//double phiscale, temp;
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.uc[0];
	//	if (phiscale < iinv.uc[cId])
	//	{
	//		phiscale = iinv.uc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.buc[cId]/(iinv.spc[cId]*phiscale);
	//		iinv.res_u += temp * temp;
	//	}

	//}
	//iinv.res_u = sqrt(iinv.res_u);

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.vc[0];
	//	if (phiscale < iinv.vc[cId])
	//	{
	//		phiscale = iinv.vc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.bvc[cId] / (iinv.spc[cId] * phiscale);
	//		iinv.res_v += temp * temp;
	//	}

	//}
	//iinv.res_v = sqrt(iinv.res_v);

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	phiscale = iinv.wc[0];
	//	if (phiscale < iinv.wc[cId])
	//	{
	//		phiscale = iinv.wc[cId];
	//	}
	//}
	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	if (iinv.spc[cId] * phiscale - 0.0 > 1e-6)
	//	{
	//		temp = iinv.bwc[cId] / (iinv.spc[cId] * phiscale);
	//		iinv.res_w += temp * temp;
	//	}

	//}
	//iinv.res_w = sqrt(iinv.res_w);


	/*for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.res_u += (iinv.buc[ug.cId]+iinv.muc[ug.cId] - iinv.ump[ug.cId]* (iinv.spu[ug.cId]))*(iinv.buc[ug.cId]+iinv.muc[ug.cId]  - iinv.ump[ug.cId] * (iinv.spu[ug.cId]));
		iinv.res_v += (iinv.bvc[ug.cId]+iinv.mvc[ug.cId] - iinv.vmp[ug.cId] * (iinv.spv[ug.cId]))*(iinv.bvc[ug.cId]+iinv.mvc[ug.cId] - iinv.vmp[ug.cId] * (iinv.spv[ug.cId]));
		iinv.res_w += (iinv.bwc[ug.cId]+iinv.mwc[ug.cId] - iinv.wmp[ug.cId] * (iinv.spw[ug.cId]))*(iinv.bwc[ug.cId]+iinv.mwc[ug.cId] - iinv.wmp[ug.cId] * (iinv.spw[ug.cId]));
	}

	iinv.res_u = sqrt(iinv.res_u);
	iinv.res_v = sqrt(iinv.res_v);
	iinv.res_w = sqrt(iinv.res_w);*/

}

void UINsInvterm::AddFlux()
{
	UnsGrid* grid = Zone::GetUnsGrid();
	MRField* res = GetFieldPointer< MRField >(grid, "res");
	int nEqu = res->GetNEqu();
	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];
		//if ( ug.lc == 0 ) std::cout << fId << std::endl;

		for (int iEqu = 0; iEqu < nEqu; ++iEqu)
		{
			(*res)[iEqu][ug.lc] -= (*iinvflux)[iEqu][ug.fId];
		}
	}

	for (int fId = ug.nBFaces; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//if ( ug.lc == 0 || ug.rc == 0 ) std::cout << fId << std::endl;

		for (int iEqu = 0; iEqu < nEqu; ++iEqu)
		{
			(*res)[iEqu][ug.lc] -= (*iinvflux)[iEqu][ug.fId];
			(*res)[iEqu][ug.rc] += (*iinvflux)[iEqu][ug.fId];
		}
  }

	//ONEFLOW::AddF2CField(res, iinvflux);
}

void UINsInvterm::CalcCorrectPresscoef()
{
	this->CalcNewMomCoe();
	for (int fId = ug.nBFaces; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CalcINsFaceCorrectPresscoef();
	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		this->CalcINsBcFaceCorrectPresscoef();
	}

	iinv.spp = 0;
	iinv.bp = 0;

	for (int fId = 0; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spp[ug.lc] += iinv.ajp[ug.fId];
		iinv.spp[ug.rc] += iinv.ajp[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];

		if (ug.fId < ug.nBFaces)
		{
			//iinv.spp[ug.rc] = 0.001;
			iinv.spp[ug.rc] = 1;
		}
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		//iinv.VdU[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spu[ug.cId] - iinv.sju[ug.cId]); //It is used to calculate the unit correction speed;
		//iinv.VdV[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spv[ug.cId] - iinv.sjv[ug.cId]);
		//iinv.VdW[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spw[ug.cId] - iinv.sjw[ug.cId]);

		iinv.VdU[ug.cId] = -(*ug.cvol)[ug.cId] / ((1+1)*iinv.spc[ug.cId]); //It is used to calculate the unit correction speed;
		iinv.VdV[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spc[ug.cId]);
		iinv.VdW[ug.cId] = -(*ug.cvol)[ug.cId] / ((1 + 1)*iinv.spc[ug.cId]);

		//iinv.spp[ug.cId] = (*ug.cvol)[ug.cId] / iinv.timestep;

		//iinv.bp[ug.cId] = iinv.bi1[ug.cId]+ iinv.bi2[ug.cId];

		int fn = (*ug.c2f)[ug.cId].size();
		if (ctrl.currTime == 0.001 && Iteration::innerSteps == 1)
		{
			iinv.sjp.resize(ug.nTCell, fn);
			iinv.sjd.resize(ug.nTCell, fn);
		}
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (ug.cId == ug.lc)
			{
				iinv.sjp[ug.cId][iFace] = -iinv.ajp[ug.fId]; //Non zero coefficient for solving pressure correction equation
				iinv.sjd[ug.cId][iFace] = ug.rc;

				//cout << "iinv.sjp=" << iinv.sjp[ug.cId][iFace] << "iinv.sjd=" << ug.rc << "\n";
			}
			else if (ug.cId == ug.rc)
			{
				iinv.sjp[ug.cId][iFace] = -iinv.ajp[ug.fId];
				iinv.sjd[ug.cId][iFace] = ug.lc;

				//cout << "iinv.sjp=" << iinv.sjp[ug.cId][iFace] << "iinv.sjd=" << ug.lc << "\n";
			}
		}
	}

	/*for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		cout << "iinv.bp=" << iinv.buc[ug.cId] << "\n";
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		cout << "iinv.spp=" << iinv.spp[ug.cId] << "\n";
	}*/

	//iinv.spp[0] = 3.996004185733362E-003;
	//iinv.spp[1] = 3.996004185733362E-003;
	//iinv.spp[2] = 3.996004185733362E-003;
	//iinv.spp[3] = 3.996004185733362E-003;
	//iinv.spp[4] = 3.996004185733362E-003;
	//iinv.spp[5] = 3.996004185733362E-003; 
	//iinv.spp[6] = 3.996004185733362E-003;
	//iinv.spp[7] = 3.996004185733362E-003;
	//iinv.spp[8] = 3.996004185733362E-003; 
	//iinv.spp[9] = 3.996004185733362E-003; 
	//iinv.spp[10] = 3.996004185733362E-003;
	//iinv.spp[11] = 3.996004185733362E-003;
	//iinv.spp[12] = 3.996004185733362E-003;
	//iinv.spp[13] = 3.996004185733362E-003; 
	//iinv.spp[14] = 3.996004185733362E-003;
	//iinv.spp[15] = 3.996003562604947E-003; 
	//	iinv.spp[16] = 3.996004185576957E-003; 
	//	iinv.spp[17] = 3.996004185733322E-003; 
	//	iinv.spp[18] = 3.996004185733362E-003; 
	//	iinv.spp[19] = 3.996004809018999E-003; 
	//	iinv.spp[20] = 3.993512901138565E-003; 
	//	iinv.spp[21] = 3.996003562135695E-003; 
	//	iinv.spp[22] = 3.996004185576840E-003;
	//	iinv.spp[23] = 3.996004809018960E-003;
	//	iinv.spp[24] = 3.998507933456209E-003;

		//for (int cId = 0; cId < 20; ++cId)
		//{
		//	ug.cId = cId;

		//	iinv.bp[ug.cId] = 0;
		//}

		//iinv.bp[20] = 9.990010315470651E-005;
		//iinv.bp[21] = 2.507471755350644E-008;
		//iinv.bp[22] = 6.309350054906251E-012;
		//iinv.bp[23] = 1.587575580783794E-015;
		//iinv.bp[24] = -9.992518418319765E-005;


}

void UINsInvterm::CalcNewMomCoe()
{
	iinv.spc = 0;

	for (int fId = 0; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.spc[ug.lc] += iinv.ai[ug.fId][0] + iinv.Fn[ug.fId];
		iinv.spc[ug.rc] += iinv.ai[ug.fId][1] + iinv.Fn[ug.fId];
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;
		iinv.spc[ug.cId] += iinv.spt[ug.cId];
	}

	//for (int cId = 0; cId < ug.nTCell; ++cId)
	//{
	//	ug.cId = cId;

	//	iinv.spu[ug.cId] = iinv.bi1[ug.cId] + iinv.bi2[ug.cId] + iinv.aku1[ug.cId] + iinv.aku2[ug.cId] + iinv.spt[ug.cId]; //The main diagonal coefficient of matrix and the principal coefficient of element of momentum equation
	//	iinv.spv[ug.cId] = iinv.bi1[ug.cId] + iinv.bi2[ug.cId] + iinv.akv1[ug.cId] + iinv.akv2[ug.cId] + iinv.spt[ug.cId];
	//	iinv.spw[ug.cId] = iinv.bi1[ug.cId] + iinv.bi2[ug.cId] + iinv.akw1[ug.cId] + iinv.akw2[ug.cId] + iinv.spt[ug.cId];
	//}

}

void UINsInvterm::CalcPressCorrectEqu()
{
	/*double rhs_p = 1e-8;
	iinv.res_p = 1;
	iinv.mp = 0;
    iinv.pp = 0;
	while (iinv.res_p >= rhs_p)
	{
		iinv.res_p = 0.0;

		for (int cId = 0; cId < ug.nCells; ++cId)
		{
			ug.cId = cId;

			iinv.ppd = iinv.pp[ug.cId];
			int fn = (*ug.c2f)[ug.cId].size();
			for (int iFace = 0; iFace < fn; ++iFace)
			{
				int fId = (*ug.c2f)[ug.cId][iFace];
				ug.fId = fId;
				if (ug.fId < ug.nBFaces) continue;

				ug.lc = (*ug.lcf)[ug.fId];
				ug.rc = (*ug.rcf)[ug.fId];
				if (ug.cId == ug.lc)
				{
					iinv.mp[ug.cId] += -iinv.sjp[ug.cId][iFace] * iinv.pp[ug.rc]; //The matrix method does not need the values of adjacent elements in Gauss Seidel iteration
				}
				else if (ug.cId == ug.rc)
				{
					iinv.mp[ug.cId] += -iinv.sjp[ug.cId][iFace] * iinv.pp[ug.lc];
				}
			}
			iinv.pp[ug.cId] = (iinv.bp[ug.cId] + iinv.mp[ug.cId]) / (iinv.spp[ug.cId]); //Pressure correction value

			iinv.res_p = MAX(iinv.res_p, abs(iinv.ppd - iinv.pp[ug.cId]));

		}

	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.pp[ug.rc] = iinv.pp[ug.lc];
	}



	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;
		(*uinsf.q)[IIDX::IIP][ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] + 0.8*iinv.pp[ug.cId];
	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
	}*/

		//Bgmres solution
	NonZero.Number = 0;

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{   
		//ug.cId = cId;                                                                  // Main unit number
		int fn = (*ug.c2f)[cId].size();                                                                 // Number of adjacent faces of element
		NonZero.Number += fn;
	}
	NonZero.Number = NonZero.Number + ug.nTCell;                                                        // Count of non-zero elements
	Rank.RANKNUMBER = ug.nTCell;                                                                        // Row and column of matrix
	Rank.COLNUMBER = 1;
	Rank.NUMBER = NonZero.Number;                                                                      // The number of nonzero elements in matrix
	Rank.Init();
	double residual_p;
	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		iinv.ppd = iinv.pp[cId];
		Rank.TempIA[0] = 0;
		int n = Rank.TempIA[cId];
		int fn = (*ug.c2f)[cId].size();
		Rank.TempIA[cId + 1] = Rank.TempIA[cId] + fn + 1;                  // The number of non-zero elements in the first n + 1 row
		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[cId][iFace];                           // Number of adjacent faces
			ug.fId = fId;
			ug.lc = (*ug.lcf)[fId];                                    // Face left unit
			ug.rc = (*ug.rcf)[fId];                                    // Face right unit
			if (cId == ug.lc)
			{
				Rank.TempA[n + iFace] = iinv.sjp[cId][iFace];          //Non diagonal element value
				Rank.TempJA[n + iFace] = ug.rc;                           //Vertical coordinates of non diagonal element
			}
			else if (cId == ug.rc)
			{
				Rank.TempA[n + iFace] = iinv.sjp[cId][iFace];          //Non diagonal element value
				Rank.TempJA[n + iFace] = ug.lc;                           //Vertical coordinates of non diagonal element
			}
		}
		Rank.TempA[n + fn] = iinv.spp[cId];                            //Main diagonal element
		Rank.TempJA[n + fn] = cId;                                        //Principal diagonal ordinate

		Rank.TempB[cId][0] = iinv.bp[cId];                             //Right end item
	}
	bgx.BGMRES();
	residual_p = Rank.residual;
	//cout << "residual_p:" << residual_p << std::endl;
	for (int cId = 0; cId < ug.nTCell; cId++)
	{
		//ug.cId = cId;
		iinv.pp[cId] = Rank.TempX[cId][0]; //Of the current momentPressure correction value
	}

	Rank.Deallocate();

	//iinv.res_p = 0;
	//iinv.res_p = MAX(iinv.res_p, abs(iinv.ppd - iinv.pp[ug.cId]));

	//boundary element
	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.pp[ug.rc] = iinv.pp[ug.lc];
	}

	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;
		(*uinsf.q)[IIDX::IIP][ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] +0.8*iinv.pp[ug.cId];
	}


	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		if (ug.bctype < 0)
		{
			false;
		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
		    (*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

		else if (ug.bctype == BC::INFLOW)
		{
			(*uinsf.q)[IIDX::IIP][ug.rc] = nscom.inflow[IIDX::IIP];
		}

		else if (ug.bctype == BC::OUTFLOW)
		{
			(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

		else if (ug.bctype == BC::POLE || ug.bctype / 10 == BC::POLE)
		{
			(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

		else if (ug.bctype == BC::EXTRAPOLATION)
		{
			(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

		else if (ug.bctype == BC::FARFIELD)
		{
			Real rin = (*uinsf.q)[IIDX::IIR][ug.lc];
			Real uin = iinv.uc[ug.lc];
			Real vin = iinv.vc[ug.lc];
			Real win = iinv.wc[ug.lc];
			Real pin = (*uinsf.q)[IIDX::IIP][ug.lc];

			gcom.xfn *= nscom.faceOuterNormal;
			gcom.yfn *= nscom.faceOuterNormal;
			gcom.zfn *= nscom.faceOuterNormal;

			Real rref = nscom.inflow[IIDX::IIR];
			Real uref = nscom.inflow[IIDX::IIU];
			Real vref = nscom.inflow[IIDX::IIV];
			Real wref = nscom.inflow[IIDX::IIW];
			Real pref = nscom.inflow[IIDX::IIP];

			Real vnref = gcom.xfn * uref + gcom.yfn * vref + gcom.zfn * wref - gcom.vfn;
			Real vnin = gcom.xfn * uin + gcom.yfn * vin + gcom.zfn * win - (*ug.vfn)[ug.fId];

			Real cref = sqrt(ABS(nscom.gama_ref * pref / rref));
			Real cin = sqrt(ABS(nscom.gama * pin / rin));

			Real gamm1 = nscom.gama - one;

			Real velin = DIST(uin, vin, win);
			//Supersonic
			if (velin > cin)
			{
				if (vnin >= 0.0)
				{
					(*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
				}
				else
				{
					(*uinsf.q)[IIDX::IIP][ug.rc] = nscom.inflow[IIDX::IIP];
				}
			}
			else
			{
				//subsonic
				Real riemp = vnin + 2.0 * cin / gamm1;
				Real riemm = vnref - 2.0 * cref / gamm1;
				Real vnb = half * (riemp + riemm);
				Real cb = fourth * (riemp - riemm) * gamm1;

				Real entr;
				if (vnb >= 0.0)
				{
					// exit
					entr = pin / pow(rin, nscom.gama);
				}
				else
				{
					entr = pref / pow(rref, nscom.gama);
				}

				Real rb = pow((cb * cb / (entr * nscom.gama)), one / gamm1);
				Real pb = cb * cb * rb / nscom.gama;

				(*uinsf.q)[IIDX::IIP][ug.rc] = pb;
			}
		}

		else
		{
		   (*uinsf.q)[IIDX::IIP][ug.rc] = (*uinsf.q)[IIDX::IIP][ug.lc];
		}

	}

	/*for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;
		iinv.pp[ug.cId] = 0;
	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.pp[ug.rc] = iinv.pp[ug.lc];
	}

for (int cId = 0; cId < ug.nCells; ++cId)
{
	ug.cId = cId;
	(*uinsf.q)[IIDX::IIP][ug.cId] = (*uinsf.q)[IIDX::IIP][ug.cId] + 0.8*iinv.pp[ug.cId];
}*/

	//for (int cId = 0; cId < ug.nTCell; cId++)
	//{
	//	iinv.pc[ug.cId] = nscom.prim[IIDX::IIP] + iinv.pp[ug.cId]; //Pressure value at the next moment
	//}
	
	/*ofstream fileres_p("residual_p.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_p << residual_p << std::endl;
	fileres_p.close();*/

}


void UINsInvterm::CalcINsPreRes()
{
	//iinv.res_p = 0;


	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		if (ug.cId == 0)
		{
			iinv.res_p = abs(iinv.bp[ug.cId]);
		}
		else
		{
			iinv.res_p = MAX(abs(iinv.bp[ug.cId]), abs(iinv.bp[ug.cId - 1]));
		}

		//iinv.res_p += (iinv.bp[ug.cId]+iinv.mp[ug.cId] - iinv.pp1[ug.cId]* (0.01+iinv.spp[ug.cId]))*(iinv.bp[ug.cId]+iinv.mp[ug.cId] - iinv.pp1[ug.cId]* (0.01+iinv.spp[ug.cId]));
	}

	//iinv.res_p = sqrt(iinv.res_p);
	//iinv.res_p = 0;
}


void UINsInvterm::UpdateFaceflux()
{
	iinv.Init();
	ug.Init();
	uinsf.Init();
	//Alloc();
	//this->CalcInvFace();  //Boundary treatment
	for (int fId = ug.nBFaces; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CalcUpdateINsFaceflux();

	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		if (fId == 10127)
		{
			int kkk = 1;
		}

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//this->PrepareFaceValue();

		this->CalcUpdateINsBcFaceflux();
	}

}

void UINsInvterm::CalcUpdateINsBcFaceflux()
{
	iinv.uuj[ug.fId] = 0; //Surface velocity correction
	iinv.vvj[ug.fId] = 0;
	iinv.wwj[ug.fId] = 0;

	iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuj[ug.fId]; //Next moment surface velocity
	iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvj[ug.fId];
	iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwj[ug.fId];

	iinv.fux = iinv.rf[ug.fId] * (gcom.xfn * iinv.uuj[ug.fId] + gcom.yfn * iinv.vvj[ug.fId] + gcom.zfn * iinv.wwj[ug.fId]) * (*ug.farea)[ug.fId];
	iinv.fq[ug.fId] = iinv.fq[ug.fId] + iinv.fux;
}


void UINsInvterm::CalcUpdateINsFaceflux()
{

	iinv.dist = (*ug.xfn)[ug.fId] * ((*ug.xcc)[ug.rc] - (*ug.xcc)[ug.lc]) + (*ug.yfn)[ug.fId] * ((*ug.ycc)[ug.rc] - (*ug.ycc)[ug.lc]) + (*ug.zfn)[ug.fId] * ((*ug.zcc)[ug.rc] - (*ug.zcc)[ug.lc]);

	iinv.uuj[ug.fId] = iinv.Vdvu[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.xfn)[ug.fId] / iinv.dist; //Surface velocity correction
	iinv.vvj[ug.fId] = iinv.Vdvv[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.yfn)[ug.fId] / iinv.dist;
	iinv.wwj[ug.fId] = iinv.Vdvw[ug.fId] * (iinv.pp[ug.lc] - iinv.pp[ug.rc]) * (*ug.zfn)[ug.fId] / iinv.dist;

	iinv.uf[ug.fId] = iinv.uf[ug.fId] + iinv.uuj[ug.fId]; //Next moment surface velocity
	iinv.vf[ug.fId] = iinv.vf[ug.fId] + iinv.vvj[ug.fId];
	iinv.wf[ug.fId] = iinv.wf[ug.fId] + iinv.wwj[ug.fId];

	iinv.fux = iinv.rf[ug.fId] * ((*ug.xfn)[ug.fId] * iinv.uuj[ug.fId] + (*ug.yfn)[ug.fId] * iinv.vvj[ug.fId] + (*ug.zfn)[ug.fId] * iinv.wwj[ug.fId]) * (*ug.farea)[ug.fId];
	iinv.fq[ug.fId] = iinv.fq[ug.fId] + iinv.fux;

}

void UINsInvterm::UpdateSpeed()
{
	this->CalcPreGrad();

	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = iinv.VdU[ug.cId] * iinv.dqqdx[ug.cId]*0.8; //Speed correction
		iinv.vv[ug.cId] = iinv.VdV[ug.cId] * iinv.dqqdy[ug.cId]*0.8;
		iinv.ww[ug.cId] = iinv.VdW[ug.cId] * iinv.dqqdz[ug.cId]*0.8;

		iinv.up[ug.cId] = iinv.uc[cId] + iinv.uu[ug.cId];  //Speed at the next moment
		iinv.vp[ug.cId] = iinv.vc[cId] + iinv.vv[ug.cId];
		iinv.wp[ug.cId] = iinv.wc[cId] + iinv.ww[ug.cId];

		(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
		(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
		(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];

	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;

		BcInfo * bcInfo = ug.bcRecord->bcInfo;

		ug.fId = bcInfo->bcFace[ug.ir][fId];
		ug.bcNameId = bcInfo->bcNameId[ug.ir][fId];

		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		if (ug.bctype < 0)
		{
			false;
		}

		else if (ug.bctype == BC::SOLID_SURFACE)
		{
			nscom.bcdtkey = 0;
			if (ug.bcNameId == -1) return; //interface
			int dd = ins_bc_data.r2d[ug.bcNameId];
			if (dd != -1)
			{
				nscom.bcdtkey = 1;
				nscom.bcflow = &ins_bc_data.dataList[dd];
			}

			if (nscom.bcdtkey == 0)
			{
				iinv.up[ug.rc] = -iinv.up[ug.lc] + 2 * gcom.vfx;
				iinv.vp[ug.rc] = -iinv.vp[ug.lc] + 2 * gcom.vfy;
				iinv.wp[ug.rc] = -iinv.wp[ug.lc] + 2 * gcom.vfz;
			}
			else
			{
				iinv.up[ug.rc] = -iinv.up[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIU];
				iinv.vp[ug.rc] = -iinv.vp[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIV];
				iinv.wp[ug.rc] = -iinv.wp[ug.lc] + 2 * (*nscom.bcflow)[IIDX::IIW];
			}
		}

		else if (ug.bctype == BC::INFLOW)
		{
			iinv.up[ug.rc] = nscom.inflow[IIDX::IIU];
			iinv.vp[ug.rc] = nscom.inflow[IIDX::IIV];
			iinv.wp[ug.rc] = nscom.inflow[IIDX::IIW];
		}

		else if (ug.bctype == BC::OUTFLOW)
		{
			iinv.up[ug.rc] = iinv.up[ug.lc];
			iinv.vp[ug.rc] = iinv.vp[ug.lc];
			iinv.wp[ug.rc] = iinv.wp[ug.lc];
		}

		else if (ug.bctype == BC::POLE || ug.bctype / 10 == BC::POLE)
		{
			iinv.up[ug.rc] = iinv.up[ug.lc];
			iinv.vp[ug.rc] = iinv.vp[ug.lc];
			iinv.wp[ug.rc] = iinv.wp[ug.lc];
		}

		else if (ug.bctype == BC::EXTRAPOLATION)
		{
			iinv.up[ug.rc] = iinv.up[ug.lc];
			iinv.vp[ug.rc] = iinv.vp[ug.lc];
			iinv.wp[ug.rc] = iinv.wp[ug.lc];
		}

		else if (ug.bctype == BC::SYMMETRY)
		{
			Real vx1 = iinv.up[ug.lc];
			Real vy1 = iinv.vp[ug.lc];
			Real vz1 = iinv.wp[ug.lc];

			Real vnRelative1 = (*ug.xfn)[ug.fId] * vx1 + (*ug.yfn)[ug.fId] * vy1 + (*ug.zfn)[ug.fId] * vz1 - (*ug.vfn)[ug.fId];

			iinv.up[ug.rc] = iinv.up[ug.lc] - two * (*ug.xfn)[ug.fId] * vnRelative1;
			iinv.vp[ug.rc] = iinv.vp[ug.lc] - two * (*ug.yfn)[ug.fId] * vnRelative1;
			iinv.wp[ug.rc] = iinv.wp[ug.lc] - two * (*ug.zfn)[ug.fId] * vnRelative1;

		}

		else if (ug.bctype == BC::SYMMETRY)
		{

		}

		else if (ug.bctype == BC::FARFIELD)
		{
			Real rin = (*uinsf.q)[IIDX::IIR][ug.lc];
			Real uin = iinv.up[ug.lc];
			Real vin = iinv.vp[ug.lc];
			Real win = iinv.wp[ug.lc];
			Real pin = (*uinsf.q)[IIDX::IIP][ug.lc];

			gcom.xfn *= nscom.faceOuterNormal;
			gcom.yfn *= nscom.faceOuterNormal;
			gcom.zfn *= nscom.faceOuterNormal;

			Real rref = nscom.inflow[IIDX::IIR];
			Real uref = nscom.inflow[IIDX::IIU];
			Real vref = nscom.inflow[IIDX::IIV];
			Real wref = nscom.inflow[IIDX::IIW];
			Real pref = nscom.inflow[IIDX::IIP];

			Real vnref = gcom.xfn * uref + gcom.yfn * vref + gcom.zfn * wref - gcom.vfn;
			Real vnin = gcom.xfn * uin + gcom.yfn * vin + gcom.zfn * win - (*ug.vfn)[ug.fId];

			Real cref = sqrt(ABS(nscom.gama_ref * pref / rref));
			Real cin = sqrt(ABS(nscom.gama * pin / rin));

			Real gamm1 = nscom.gama - one;

			Real velin = DIST(uin, vin, win);

			//Supersonic
			if (velin > cin)
			{
				if (vnin >= 0.0)
				{
					iinv.up[ug.rc] = iinv.up[ug.lc];
					iinv.vp[ug.rc] = iinv.vp[ug.lc];
					iinv.wp[ug.rc] = iinv.wp[ug.lc];
				}
				else
				{
					iinv.up[ug.rc] = nscom.inflow[IIDX::IIU];
					iinv.vp[ug.rc] = nscom.inflow[IIDX::IIV];
					iinv.wp[ug.rc] = nscom.inflow[IIDX::IIW];
				}
			}
			else
			{
				//subsonic
				Real riemp = vnin + 2.0 * cin / gamm1;
				Real riemm = vnref - 2.0 * cref / gamm1;
				Real vnb = half * (riemp + riemm);
				Real cb = fourth * (riemp - riemm) * gamm1;

				Real vtx, vty, vtz, entr;
				if (vnb >= 0.0)
				{
					// exit
					entr = pin / pow(rin, nscom.gama);

					vtx = uin - gcom.xfn * vnin;
					vty = vin - gcom.yfn * vnin;
					vtz = win - gcom.zfn * vnin;
				}
				else
				{
					//inlet
					entr = pref / pow(rref, nscom.gama);
					vtx = uref - gcom.xfn * vnref;
					vty = vref - gcom.yfn * vnref;
					vtz = wref - gcom.zfn * vnref;
				}

				Real rb = pow((cb * cb / (entr * nscom.gama)), one / gamm1);
				Real ub = vtx + gcom.xfn * vnb;
				Real vb = vty + gcom.yfn * vnb;
				Real wb = vtz + gcom.zfn * vnb;
				Real pb = cb * cb * rb / nscom.gama;


				iinv.up[ug.rc] = ub;
				iinv.vp[ug.rc] = vb;
				iinv.wp[ug.rc] = wb;

			}
		}

		else if (ug.bctype == BC::OVERSET)
		{
			;
		}

		else if (ug.bctype == BC::GENERIC_2)
		{

			;

		}


		(*uinsf.q)[IIDX::IIU][ug.rc] = iinv.up[ug.rc];
		(*uinsf.q)[IIDX::IIV][ug.rc] = iinv.vp[ug.rc];
		(*uinsf.q)[IIDX::IIW][ug.rc] = iinv.wp[ug.rc];

	}



	/*for (int cId = ug.nCells; cId < ug.nTCell; ++cId)
	{
		ug.cId = cId;

		iinv.uu[ug.cId] = 0; //Speed correction
		iinv.vv[ug.cId] = 0;
		iinv.ww[ug.cId] = 0;

		iinv.up[ug.cId] = iinv.uc[cId] + iinv.uu[ug.cId];  //Speed at the next moment
		iinv.vp[ug.cId] = iinv.vc[cId] + iinv.vv[ug.cId];
		iinv.wp[ug.cId] = iinv.wc[cId] + iinv.ww[ug.cId];

		(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
		(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
		(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];
	}*/
}

void UINsInvterm::UpdateINsRes()
{
	/*iinv.remax_V = 0;
	iinv.remax_pp = 0;

	for (int fId = 0; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
	}

	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;
		iinv.res_V[ug.cId] = 10*iinv.bp[ug.cId];

		iinv.remax_V = MAX(iinv.remax_V, abs(iinv.res_V[ug.cId]));
		iinv.remax_pp = MAX(iinv.remax_pp, abs(iinv.pp[ug.cId]));

	}
	cout << "iinv.remax_V:" << iinv.remax_V << std::endl;
	cout << "iinv.remax_pp:" << iinv.remax_pp << std::endl;
	cout <<"innerSteps:"<< Iteration::innerSteps<< std::endl;
	//cout << "outerSteps:" << Iteration::outerSteps << std::endl;

	ofstream fileres_vv("residual_vv.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_vv << iinv.remax_V << std::endl;
	fileres_vv.close();
	

	ofstream fileres_pp("residual_pp.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_pp << iinv.remax_pp << std::endl;
	fileres_pp.close();*/






	iinv.remax_up = 0;
	iinv.remax_vp = 0;
	iinv.remax_wp = 0;
	iinv.remax_pp = 0;

	for (int fId = 0; fId < ug.nFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		iinv.bp[ug.lc] += -iinv.fq[ug.fId];
		iinv.bp[ug.rc] += iinv.fq[ug.fId];
	}

	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;

		int fn = (*ug.c2f)[ug.cId].size();

		for (int iFace = 0; iFace < fn; ++iFace)
		{
			int fId = (*ug.c2f)[ug.cId][iFace];
			ug.fId = fId;
			ug.lc = (*ug.lcf)[ug.fId];
			ug.rc = (*ug.rcf)[ug.fId];

			if (ug.cId == ug.lc)
			{
				iinv.mu[ug.cId] += -iinv.ai[ug.fId][0] * (iinv.up[ug.rc] - iinv.uc[ug.rc]);  //The flux of the element surface adjacent to the main element in the momentum equation
				iinv.mv[ug.cId] += -iinv.ai[ug.fId][0] * (iinv.vp[ug.rc] - iinv.vc[ug.rc]);
				iinv.mw[ug.cId] += -iinv.ai[ug.fId][0] * (iinv.wp[ug.rc] - iinv.wc[ug.rc]);
				//iinv.mpp[ug.cId] += -iinv.ajp[ug.fId] * iinv.pp[ug.rc];
			}
			else if (ug.cId == ug.rc)
			{
				iinv.mu[ug.cId] += -iinv.ai[ug.fId][1] * (iinv.up[ug.lc] - iinv.uc[ug.lc]);  //The flux of the element surface adjacent to the main element in the momentum equation
				iinv.mv[ug.cId] += -iinv.ai[ug.fId][1] * (iinv.vp[ug.lc] - iinv.vc[ug.lc]);
				iinv.mw[ug.cId] += -iinv.ai[ug.fId][1] * (iinv.wp[ug.lc] - iinv.wc[ug.lc]);
				//iinv.mpp[ug.cId] += -iinv.ajp[ug.fId] * iinv.pp[ug.lc];
			}
		}

		iinv.mua[ug.cId] = iinv.spc[ug.cId] * (iinv.up[ug.cId] - iinv.uc[ug.cId]) + iinv.mu[ug.cId];
		iinv.mva[ug.cId] = iinv.spc[ug.cId] * (iinv.vp[ug.cId] - iinv.vc[ug.cId]) + iinv.mv[ug.cId];
		iinv.mwa[ug.cId] = iinv.spc[ug.cId] * (iinv.wp[ug.cId] - iinv.wc[ug.cId]) + iinv.mw[ug.cId];
		//iinv.mppa[ug.cId] = iinv.spp[ug.cId] * iinv.pp[ug.cId] + iinv.mpp[ug.cId];

		iinv.res_up[ug.cId] = iinv.mua[ug.cId] * iinv.mua[ug.cId];
		iinv.res_vp[ug.cId] = iinv.mva[ug.cId] * iinv.mva[ug.cId];
		iinv.res_wp[ug.cId] = iinv.mwa[ug.cId] * iinv.mwa[ug.cId];
		iinv.res_pp[ug.cId] = iinv.bp[ug.cId] * iinv.bp[ug.cId];

		iinv.remax_up += iinv.res_up[ug.cId];
		iinv.remax_vp += iinv.res_vp[ug.cId];
		iinv.remax_wp += iinv.res_wp[ug.cId];
		iinv.remax_pp += iinv.res_pp[ug.cId];
	}


	iinv.remax_up = sqrt(iinv.remax_up);
	iinv.remax_vp = sqrt(iinv.remax_vp);
	iinv.remax_wp = sqrt(iinv.remax_wp);
	iinv.remax_pp = sqrt(iinv.remax_pp);

	cout << "iinv.remax_up:" << iinv.remax_up << std::endl;
	cout << "iinv.remax_vp:" << iinv.remax_vp << std::endl;
	cout << "iinv.remax_wp:" << iinv.remax_wp << std::endl;
	cout << "iinv.remax_pp:" << iinv.remax_pp << std::endl;


	ofstream fileres_up("residual_up.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_up << iinv.remax_up << std::endl;
	fileres_up.close();


	ofstream fileres_vp("residual_vp.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_vp << iinv.remax_vp << std::endl;
	fileres_vp.close();

	ofstream fileres_wp("residual_wp.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_wp << iinv.remax_wp << std::endl;
	fileres_wp.close();

	ofstream fileres_pp("residual_pp.txt", std::ios::app);
	//fileres_p << "residual_p:" <<residual_p << std::endl;
	fileres_pp << iinv.remax_pp << std::endl;
	fileres_pp.close();


}


void UINsInvterm::CalcPreGrad()
{
	iinv.dqqdx = 0;
	iinv.dqqdy = 0;
	iinv.dqqdz = 0;

	for (int fId = 0; fId < ug.nFaces; ++fId)
	{
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
		Real delta = 1.0 / (delt1 + delt2 + SMALL);

		Real cl = delt2 * delta;
		Real cr = delt1 * delta;
		//if (ug.fId < ug.nBFaces)
		//{
		//	iinv.value[ug.fId] = iinv.pp[ug.lc] + iinv.pp[ug.rc];
		//}
		//else
		//{
		iinv.value = cl * iinv.pp[ug.lc] + cr * iinv.pp[ug.rc];
		//}

		Real fnxa = (*ug.xfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnya = (*ug.yfn)[ug.fId] * (*ug.farea)[ug.fId];
		Real fnza = (*ug.zfn)[ug.fId] * (*ug.farea)[ug.fId];

		iinv.dqqdx[ug.lc] += fnxa * iinv.value;
		iinv.dqqdy[ug.lc] += fnya * iinv.value;
		iinv.dqqdz[ug.lc] += fnza * iinv.value;

		if (ug.fId < ug.nBFaces) continue;

		iinv.dqqdx[ug.rc] += -fnxa * iinv.value;
		iinv.dqqdy[ug.rc] += -fnya * iinv.value;
		iinv.dqqdz[ug.rc] += -fnza * iinv.value;
	}

	for (int cId = 0; cId < ug.nCells; ++cId)
	{
		ug.cId = cId;
		Real ovol = one / (*ug.cvol)[ug.cId];
		iinv.dqqdx[ug.cId] *= ovol;
		iinv.dqqdy[ug.cId] *= ovol;
		iinv.dqqdz[ug.cId] *= ovol;
	}

	for (int fId = 0; fId < ug.nBFaces; ++fId)
	{
		ug.fId = fId;
		ug.lc = (*ug.lcf)[ug.fId];
		ug.rc = (*ug.rcf)[ug.fId];

		//if (ug.rc > ug.nCells)
		//{
		iinv.dqqdx[ug.rc] = iinv.dqqdx[ug.lc];
		iinv.dqqdy[ug.rc] = iinv.dqqdy[ug.lc];
		iinv.dqqdz[ug.rc] = iinv.dqqdz[ug.lc];
		//}

	}

}


void UINsInvterm::Alloc()
{
	//iinvflux = new MRField(nscom.nEqu, ug.nFaces);
}

void UINsInvterm::DeAlloc()
{
	//delete iinvflux;
}


void UINsInvterm::ReadTmp()
{
	static int iii = 0;
	if (iii) return;
	iii = 1;
	fstream file;
	file.open("nsflow.dat", std::ios_base::in | std::ios_base::binary);
	if (!file) exit(0);

	uinsf.Init();

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		for (int iEqu = 0; iEqu < 5; ++iEqu)
		{
			file.read(reinterpret_cast<char*>(&(*uinsf.q)[iEqu][cId]), sizeof(double));
		}
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.visl)[0][cId]), sizeof(double));
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.vist)[0][cId]), sizeof(double));
	}

	vector< Real > tmp1(ug.nTCell), tmp2(ug.nTCell);

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		tmp1[cId] = (*uinsf.timestep)[0][cId];
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		file.read(reinterpret_cast<char*>(&(*uinsf.timestep)[0][cId]), sizeof(double));
	}

	for (int cId = 0; cId < ug.nTCell; ++cId)
	{
		tmp2[cId] = (*uinsf.timestep)[0][cId];
	}

	turbcom.Init();
	uturbf.Init();
	for (int iCell = 0; iCell < ug.nTCell; ++iCell)
	{
		for (int iEqu = 0; iEqu < turbcom.nEqu; ++iEqu)
		{
			file.read(reinterpret_cast<char*>(&(*uturbf.q)[iEqu][iCell]), sizeof(double));
		}
	}
	file.close();
	file.clear();
}



EndNameSpace
