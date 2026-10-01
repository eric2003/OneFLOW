/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include <memory>


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
	auto uINsBcSolver = std::make_unique<UINsBcSolver>();
	uINsBcSolver->Init();
	uINsBcSolver->CalcBc();
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
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->CalcINsTimeStep();
}

void INsPreflux()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->CalcINsPreflux();
}

void INsCalcInv()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->CalcInvcoff();
}

void INsCalcVis()
{
	auto uINsVisterm = std::make_unique<UINsVisterm>();
	uINsVisterm->CalcViscoff();
}

void INsCalcUnstead()
{
	auto uINsVisterm = std::make_unique<UINsVisterm>();
	uINsVisterm->CalcUnsteadcoff();
}

void INsCalcSrc()
{
	auto uINsVisterm = std::make_unique<UINsVisterm>();
	uINsVisterm->CalcINsSrc();
}

void INsMomPre()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->MomPre();
}

void INsCalcFaceflux()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->CalcFaceflux();
}

void INsCorrectPresscoef()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->CalcCorrectPresscoef();
}

void INsCalcPressCorrectEquandUpdatePress()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->CalcPressCorrectEqu();
}

void INsUpdateFaceflux()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->UpdateFaceflux();
}

void INsCalcSpeedCorrectandUpdateSpeed()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->UpdateSpeed();
}

void INsUpdateRes()
{
	auto uINsInvterm = std::make_unique<UINsInvterm>();
	uINsInvterm->UpdateINsRes();
}

//void INsCorrectSpeed()
//{
//	auto uINsInvterm = std::make_unique<UINsInvterm>();
//	uINsInvterm->CalcCorrectSpeed();
//
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
//	auto uinsUnsteady = std::make_unique<UINsUnsteady>();
//	uinsUnsteady->CalcDualTimeSrc();
//
//}


EndNameSpace
