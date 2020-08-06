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

#include "NsRhs.h"
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
#include "UNsInvFlux.h"
#include "UNsVisFlux.h"
#include "UNsUnsteady.h"
#include "Ctrl.h"

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

NsRhs::NsRhs()
{
    ;
}

NsRhs::~NsRhs()
{
    ;
}

void NsRhs::UpdateResiduals()
{
	NsCalcRHS();
}

void NsCalcBc()
{
	UNsBcSolver * uNsBcSolver = new UNsBcSolver();
	uNsBcSolver->Init();
	uNsBcSolver->CalcBc();
	delete uNsBcSolver;
}

void NsCalcGamaT(int flag)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	ug.Init();
	unsf.Init();
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
			Real & density = (*unsf.q)[IDX::IR][cId];
			Real & pressure = (*unsf.q)[IDX::IP][cId];
			(*unsf.gama)[0][cId] = nscom.gama_ref;
			(*unsf.tempr)[IDX::ITT][cId] = pressure / (nscom.statecoef * density * oamw);
		}
	}
}

void NsCalcRHS()
{
	NsCalcInvFlux();

	NsCalcVisFlux();

	NsCalcSrcFlux();
}


void NsCalcInvFlux()
{
	UNsInvFlux * uNsInvFlux = new UNsInvFlux();
	uNsInvFlux->CalcFlux();
	delete uNsInvFlux;
}

void NsCalcVisFlux()
{
	UNsVisFlux * uNsVisFlux = new UNsVisFlux();
	uNsVisFlux->CalcFlux();
	delete uNsVisFlux;
}

void NsCalcSrcFlux()
{
	if (nscom.chemModel == 1)
	{
		NsCalcChemSrc();
	}

	NsCalcTurbEnergy();

	//dual time step source
	if ( ctrl.idualtime == 1 )
	{
		NsCalcDualTimeStepSrc();
	}
}

void NsCalcChemSrc()
{
	;
}

void NsCalcTurbEnergy()
{
	;
}

void NsCalcDualTimeStepSrc()
{
	UNsUnsteady * unsUnsteady = new UNsUnsteady();
	unsUnsteady->CalcDualTimeSrc();
	delete unsUnsteady;
}

EndNameSpace