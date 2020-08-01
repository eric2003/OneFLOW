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

#include "INsSolverImp.h"
#include "SolverImp.h"
#include "UITimestep.h"
#include "Zone.h"
#include "GridState.h"
#include "UINsVisterm.h"
#include "Rhs.h"
#include "UCom.h"
#include "CmxTask.h"
#include "Iteration.h"
#include "SolverDef.h"
#include "Ctrl.h"
#include "LaminarPlate.h"
#include "TurbPlate.h"
#include "SolverRegister.h"
#include "DataBase.h"

BeginNameSpace( ONEFLOW )

void RegisterINsFunc()
{
    REGISTER_DATA_CLASS( INsInitFinal );
    REGISTER_DATA_CLASS( INsVisual );
    REGISTER_DATA_CLASS( INsCmpTimeStep );
    REGISTER_DATA_CLASS( INsUpdateResiduals );
    REGISTER_DATA_CLASS( INsImplicitMethod );
    REGISTER_DATA_CLASS( INsPostprocess );
    REGISTER_DATA_CLASS( INsFinalPostprocess );
    REGISTER_DATA_CLASS( INsInitSolver );
    REGISTER_DATA_CLASS( INsCmpBoundary );
    REGISTER_DATA_CLASS( IDumpHeatFluxCoeff );

	REGISTER_DATA_CLASS( INsCmpTurb );
	REGISTER_DATA_CLASS(INsCmpHeat);
}

void INsInitFinal( StringField & data )
{
    INSCmpGamaT( F_INNER );
    //ICmpLaminarViscosity( F_INNER );
    INsCmpBc();
    INSCmpGamaT( F_GHOST );
    //ICmpLaminarViscosity( F_GHOST );

    Grid * grid = Zone::GetGrid();

    if ( Zone::GetCGrid( grid ) )
    {
        //RestrictAllQ( NS_SOLVER, FLOW_FIELD_INDEX );

        GridState::gridLevel += 1;

		NsCalcBc;

        GridState::gridLevel -= 1;
    }
}

void INsVisual( StringField & data )
{
    ;
}

void INsCmpBoundary( StringField & data )
{
    INSCmpGamaT( F_INNER );
   // ICmpLaminarViscosity( F_INNER );
    INsCmpBc();
    INSCmpGamaT( F_GHOST );
   // ICmpLaminarViscosity( F_GHOST );
}

void INsCmpTimeStep( StringField & data )
{
    UITimestep * uTimestep = new UITimestep();
    uTimestep->CmpTimestep();
    delete uTimestep;
}

void INsUpdateResiduals( StringField & data )
{
    Rhs * rhs = new Rhs();
    rhs->UpdateINsResiduals();
    delete rhs;
}

void INsImplicitMethod( StringField & data )
{

    ;
}

void INsPostprocess( StringField & data )
{
    //After every Iter, the first thing to consider is communication.
    CommInterfaceData();

    //The solution and output of residuals need to be judged logically.
    if ( Iteration::ResOk() )
    {
        ONEFLOW::AddCmdToList( "DUMP_RESIDUAL" );
    }

    //The solution and output of aerodynamic force need to be judged logically.

    if ( Iteration::ForceOk() )
    {
        ONEFLOW::AddCmdToList( "DUMP_AERODYNAMIC" );
    }

	int startStrategy = ONEFLOW::GetDataValue< int >("startStrategy");
	if (startStrategy == 2)
	{
		;
	}
	else
	{
		if (!Iteration::InnerOk()) return;

		ONEFLOW::AddCmdToList("UPDATE_UNSTEADY_FLOW");
	}

	
	if (startStrategy == 2)
	{
		if (Iteration::innerSteps % Iteration::nVisualSave == 0)
		{
			ONEFLOW::AddCmdToList("VISUALIZATION");
			ONEFLOW::AddCmdToList("DUMP_PRESSURE_COEFF");
			ONEFLOW::AddCmdToList("DUMP_HEATFLUX_COEFF");
		}
	}
	else
	{
		if (Iteration::outerSteps % Iteration::nVisualSave == 0)
		{
			ONEFLOW::AddCmdToList("VISUALIZATION");
			ONEFLOW::AddCmdToList("DUMP_PRESSURE_COEFF");
			ONEFLOW::AddCmdToList("DUMP_HEATFLUX_COEFF");
			//if ( ctrl.idump == 1 )
			//{
			//    ONEFLOW::AddCmdToList( "DUMP_LAMINAR_PLATE" );
			//}
			//else if ( ctrl.idump == 2 )
			//{
			//    ONEFLOW::AddCmdToList( "DUMP_TURB_PLATE" );
			//}
		}
	}

	if (startStrategy == 2)
	{
		if (Iteration::innerSteps % Iteration::nFieldSave == 0)
		{
			ONEFLOW::AddCmdToList("DUMP_RESTART");
		}
	}
	else
	{
		if (Iteration::outerSteps % Iteration::nFieldSave == 0)
		{
			ONEFLOW::AddCmdToList("DUMP_RESTART");
		}
	}

}

void INsFinalPostprocess( StringField & data )
{
    ONEFLOW::AddCmdToList( "DUMP_RESTART"        );
    ONEFLOW::AddCmdToList( "DUMP_AERODYNAMIC"    );
    ONEFLOW::AddCmdToList( "DUMP_PRESSURE_COEFF" );
    ONEFLOW::AddCmdToList( "DUMP_HEATFLUX_COEFF" );
    ONEFLOW::AddCmdToList( "VISUALIZATION"       );
}

void INsInitSolver( StringField & data )
{
    ONEFLOW::CommInterfaceData();
}

void IDumpHeatFluxCoeff( StringField & data )
{
    ;
}

void INsCmpTurb(StringField & data)
{
	;
}

void INsCmpHeat(StringField & data)
{
	;
}

SolverRegData insReg;
SolverRegData * GetINsReg()
{
    insReg.sTid = INC_NS_SOLVER;
    insReg.func = & RegisterINsFunc;
    insReg.solverName = "ins";
    insReg.baseKind = SOLVER_BASED;
    insReg.dataFlag = WITH_DATA;
    return & insReg;
}
  REGISTER_REG_DATA( GetINsReg );

EndNameSpace