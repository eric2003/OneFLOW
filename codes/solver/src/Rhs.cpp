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

#include "Rhs.h"
#include "UNsBcSolver.h"
#include "Zone.h"
#include "DataBase.h"
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

Rhs::Rhs()
{
    ;
}

Rhs::~Rhs()
{
    ;
}

void Rhs::UpdateNsResiduals()
{
    NsCompRHS();
}

void NsCompBc()
{
    UNsBcSolver * uNsBcSolver = new UNsBcSolver();
    uNsBcSolver->Init();
    uNsBcSolver->CalcBc();
    delete uNsBcSolver;
}

void NsCompGamaT( int flag )
{
    UnsGrid * grid = Zone::GetUnsGrid();

    ug.Init();
    unsf.Init();
    ug.SetStEd( flag );

    if ( nscom.chemModel == 1 )
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
        for ( int cId = ug.ist; cId < ug.ied; ++ cId )
        {
            Real & density  = ( * unsf.q )[ IDX::IR ][ cId ];
            Real & pressure = ( * unsf.q )[ IDX::IP ][ cId ];
            ( * unsf.gama )[ 0 ][ cId ] = nscom.gama_ref;
            ( * unsf.tempr )[ IDX::ITT ][ cId ] = pressure / ( nscom.statecoef * density * oamw );
        }
    }
}

void NsCompRHS()
{
    NsCompInvFlux();

    NsCompVisFlux();

    NsCompSrcFlux();
}


void NsCompInvFlux()
{
    UNsInvFlux * uNsInvFlux = new UNsInvFlux();
    uNsInvFlux->CalcFlux();
    delete uNsInvFlux;
}

void NsCompVisFlux()
{
    UNsVisFlux * uNsVisFlux = new UNsVisFlux();
    uNsVisFlux->CalcFlux();
    delete uNsVisFlux;
}

void NsCompSrcFlux()
{
    if ( nscom.chemModel == 1 )
    {
        NsCompChemSrc();
    }

    NsCompTurbEnergy();

    //dual time step source
    if ( ctrl.idualtime == 1 )
    {
        NsCompDualTimeStepSrc();
    }
}

void NsCompChemSrc()
{
    ;
}

void NsCompTurbEnergy()
{
    ;
}

void NsCompDualTimeStepSrc()
{
    UNsUnsteady * unsUnsteady = new UNsUnsteady();
    unsUnsteady->CmpDualTimeSrc();
    delete unsUnsteady;
}

void Rhs::UpdateINsResiduals()
{
	INsCompRHS();
}

void INsCompBc()
{
	UINsBcSolver * uINsBcSolver = new UINsBcSolver();
	uINsBcSolver->Init();
	uINsBcSolver->CalcBc();
	delete uINsBcSolver;
}

void INsCompGamaT(int flag)
{
	UnsGrid * grid = Zone::GetUnsGrid();

	ug.Init();
	uinsf.Init();
	ug.SetStEd(flag);

	if (inscom.chemModel == 1)
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
			( * uinsf.gama )[ 0 ][ cId ] = inscom.gama_ref;
			//( * uinsf.tempr )[ IIDX::IITT ][ cId ] = pressure / ( inscom.statecoef * density * oamw );
		}
	}
}

void INsCompRHS()
{
	INsCompInv();

	INsCompVis();

	INsCompSrc();

	//INsMomPred();  //需要解动量方程组

	INsCompFaceflux();

	INsCorrectPresscoef();

	INsCompPressCorrectEquandUpdatePress();  //需要解压力修正方程组，增设单元修正压力未知量

	INsCompSpeedCorrectandUpdateSpeed();  //需要先增设界面修正速度未知量并进行求解,更新单元速度和压力

	INsUpdateFaceflux();   //更新界面流量

}

void INsCompInv()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpInvcoff();
	delete uINsInvterm;
}

void INsCompVis()
{
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CalcViscoff();
	delete uINsVisterm;
}

void INsCompSrc()
{
	//if (inscom.chemModel == 1)
	//{
	//	INsCompChemSrc();
	//}

	//INsCompTurbEnergy();

	//dual time step source
	//if (ctrl.idualtime == 1)
	//{
	//	INsCompDualTimeStepSrc();
	//}
	UINsVisterm * uINsVisterm = new UINsVisterm();
	uINsVisterm->CalcSrc();
	delete uINsVisterm;
}

void INsMomPred()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->MomPred();
	delete uINsInvterm;
}

void INsCompFaceflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpFaceflux();
	delete uINsInvterm;
}

void INsCorrectPresscoef()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpCorrectPresscoef();
	delete uINsInvterm;
}

void INsCompPressCorrectEquandUpdatePress()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->CmpPressCorrectEqu();
	delete uINsInvterm;
}

void INsUpdateFaceflux()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateFaceflux();
	delete uINsInvterm;
}

void INsCompSpeedCorrectandUpdateSpeed()
{
	UINsInvterm * uINsInvterm = new UINsInvterm();
	uINsInvterm->UpdateSpeed();
	delete uINsInvterm;
}

//void INsCorrectSpeed()
//{
//	UINsInvterm * uINsInvterm = new UINsInvterm();
//	uINsInvterm->CmpCorrectSpeed();
//	delete uINsInvterm;
//}





void INsCompChemSrc()
{
	;
}

void INsCompTurbEnergy()
{
	;
}

//void INsCompDualTimeStepSrc()
//{
//	UINsUnsteady * uinsUnsteady = new UINsUnsteady();
//	uinsUnsteady->CmpDualTimeSrc();
//	delete uinsUnsteady;
//}


EndNameSpace