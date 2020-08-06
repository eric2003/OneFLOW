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

#include "INsRhs.h"
#include "UNsBcSolver.h"
#include "Zone.h"
#include "DataBase.h"
#include "Iteration.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "UnsGrid.h"
#include "NsCom.h"
#include "NsIdx.h"
#include "Ctrl.h"

//#include "UINsCorrectPress.h"
//#include "UINsCorrectSpeed.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "INsCom.h"
#include "INsIdx.h"
#include "UINsInvterm.h"
#include "UINsVisterm.h"
//#include "UINsUnsteady.h"
#include "UINsBcSolver.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

INsRhs::INsRhs()
{
    ;
}

INsRhs::~INsRhs()
{
    ;
}

void INsRhs::UpdateResiduals()
{
	INsCalcRHS();
}

void INsCalcBc()
{
	UINsBcSolver * uINsBcSolver = new UINsBcSolver();
	uINsBcSolver->Init();
	uINsBcSolver->CalcBc();
	delete uINsBcSolver;
}

void INsCalcGamaT(int flag)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	ug.Init();
	uinsf.Init();
	ug.SetStEd(flag);

	if (nscom.chemModel == 1)
	{
	}
	else
	{
		//if ( ZoneState::zid == 0 )
		//{
		//    cout << " ug.ist = " << ug.ist  << " ug.ied = " << ug.ied << "\n";
		//    int kkk = 1;
		//}
		Real oamw = one;
		for (int cId = ug.ist; cId < ug.ied; ++cId)
		{
			Real & density = ( * uinsf.q )[ IIDX::IIR ][ cId ];
			Real & pressure = ( * uinsf.q )[ IIDX::IIP ][ cId ];

			( * uinsf.gama )[ 0 ][ cId ] = nscom.gama_ref;
			//( * uinsf.tempr )[ IIDX::IITT ][ cId ] = pressure / ( nscom.statecoef * density * oamw );
			(*uinsf.tempr)[IIDX::IITT][cId] = 0;
		}
	}
}

void INsCalcRHS()
{
	INsCalcTimeStep();

	INsPreflux();

	INsCalcInv(); //计算对流项

	INsCalcVis(); //计算扩散项

	INsCalcUnstead(); //计算非稳态项

	INsCalcSrc(); //计算源项和动量方程系数

	INsMomPre(); //求解动量方程

	INsCalcFaceflux(); //计算界面流量

	INsCorrectPresscoef(); //计算压力修正方程系数

	INsCalcPressCorrectEquandUpdatePress();  //需要解压力修正方程组，增设单元修正压力未知量

	INsCalcSpeedCorrectandUpdateSpeed();  //需要先增设界面修正速度未知量并进行求解,更新单元速度和压力

	INsUpdateFaceflux();   //更新界面流量

	INsUpdateRes();
}

void INsCalcTimeStep()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CalcINsTimeStep();
	delete uINsInvterm;
}

void INsPreflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CalcINsPreflux();
	delete uINsInvterm;
}

void INsCalcInv()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CalcInvcoff();
	delete uINsInvterm;
}

void INsCalcVis()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CalcViscoff();
	delete uINsVisterm;
}

void INsCalcUnstead()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CalcUnsteadcoff();
	delete uINsVisterm;
}

void INsCalcSrc()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CalcINsSrc();
	delete uINsVisterm;
}

void INsMomPre()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->MomPre();
	delete uINsInvterm;
}

void INsCalcFaceflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CalcFaceflux();
	delete uINsInvterm;
}

void INsCorrectPresscoef()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CalcCorrectPresscoef();
	delete uINsInvterm;
}

void INsCalcPressCorrectEquandUpdatePress()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CalcPressCorrectEqu();
	delete uINsInvterm;
}

void INsUpdateFaceflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateFaceflux();
	delete uINsInvterm;
}

void INsCalcSpeedCorrectandUpdateSpeed()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateSpeed();
	delete uINsInvterm;
}

void INsUpdateRes()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateINsRes();
	delete uINsInvterm;
}

//void INsCorrectSpeed()
//{
//	UINsInvterm * uINsInvterm = new UINsInvterm();
//	uINsInvterm->CalcCorrectSpeed();
//	delete uINsInvterm;
//}



void INsCalcChemSrc()
{
	;
}

void INsCalcTurbEnergy()
{
	;
}

//void INsCalcDualTimeStepSrc()
//{
//	UINsUnsteady * uinsUnsteady = new UINsUnsteady();
//	uinsUnsteady->CalcDualTimeSrc();
//	delete uinsUnsteady;
//}


EndNameSpace