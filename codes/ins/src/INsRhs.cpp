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
		//    std::cout << " ug.ist = " << ug.ist  << " ug.ied = " << ug.ied << "\n";
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

	INsCalcInv(); //Calculating the convective term

	INsCalcVis(); //Calculate diffusion term

	INsCalcUnstead(); //Calculating the unsteady term

	INsCalcSrc(); //Calculate the source term and momentum equation coefficients

	INsMomPre(); //Solving momentum equation

	INsCalcFaceflux(); //Calculation of interface flow

	INsCorrectPresscoef(); //Calculate the coefficient of pressure correction equation

	INsCalcPressCorrectEquandUpdatePress();  //It is necessary to solve the pressure correction equations and add a new element to correct the unknown pressure

	INsCalcSpeedCorrectandUpdateSpeed();  //It is necessary to add the interface to correct the unknown velocity and solve the problem, and update the unit velocity and pressure

	INsUpdateFaceflux();   //Update interface Flux

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
