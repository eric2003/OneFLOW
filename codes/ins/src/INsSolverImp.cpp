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

#include "INsSolverImp.h"
#include "SolverImp.h"
#include "UTimeStep.h"
#include "Zone.h"
#include "GridState.h"
#include "UINsVisterm.h"
#include "INsRhs.h"
#include "UCom.h"
#include "CmxTask.h"
#include "CmxTaskNames.h"
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
    REGISTER_DATA_CLASS( INsCalcTimeStep );
    REGISTER_DATA_CLASS( INsUpdateResiduals );
    REGISTER_DATA_CLASS( INsImplicitMethod );
    REGISTER_DATA_CLASS( INsPostprocess );
    REGISTER_DATA_CLASS( INsFinalPostprocess );
    REGISTER_DATA_CLASS( INsInitSolver );
    REGISTER_DATA_CLASS( INsCalcBoundary );
    REGISTER_DATA_CLASS( IDumpHeatFluxCoeff );

	REGISTER_DATA_CLASS( INsCalcTurb );
	REGISTER_DATA_CLASS( INsCalcHeat );
}

void INsInitFinal( StringField & data )
{
    INsCalcGamaT( F_INNER );
    //ICalcLaminarViscosity( F_INNER );
    INsCalcBc();
    INsCalcGamaT( F_GHOST );
    //ICalcLaminarViscosity( F_GHOST );

    Grid * grid = Zone::GetGrid();

    if ( Zone::GetCGrid( grid ) )
    {
        //RestrictAllQ( NS_SOLVER, FLOW_FIELD_INDEX );

        GridState::gridLevel += 1;

		INsCalcBc();

        GridState::gridLevel -= 1;
    }
}

void INsVisual( StringField & data )
{
    ;
}

void INsCalcBoundary( StringField & data )
{
    INsCalcGamaT( F_INNER );
   // ICalcLaminarViscosity( F_INNER );
    INsCalcBc();
    INsCalcGamaT( F_GHOST );
   // ICalcLaminarViscosity( F_GHOST );
}

void INsCalcTimeStep( StringField & data )
{
    UTimeStep * uTimeStep = new UTimeStep();
    uTimeStep->CalcTimeStep();
    delete uTimeStep;
}

void INsUpdateResiduals( StringField & data )
{
    Rhs * rhs = new INsRhs();
    rhs->UpdateResiduals();
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
        ONEFLOW::AddCmdToList( kDumpResidualTaskName );
    }

    //The solution and output of aerodynamic force need to be judged logically.

    if ( Iteration::ForceOk() )
    {
        ONEFLOW::AddCmdToList( kDumpAerodynamicTaskName );
    }

    int startStrategy = ctrl.startStrategy;
	if (startStrategy == 2)
	{
		;
	}
	else
	{
		if (!Iteration::InnerOk()) return;

		ONEFLOW::AddCmdToList( kUpdateUnsteadyFlowTaskName );
	}

	
	if (startStrategy == 2)
	{
		if (Iteration::innerSteps % Iteration::nVisualSave == 0)
		{
			ONEFLOW::AddCmdToList( kVisualizationTaskName     );
			ONEFLOW::AddCmdToList( kDumpPressureCoeffTaskName );
			ONEFLOW::AddCmdToList( kDumpHeatfluxCoeffTaskName );
		}
	}
	else
	{
		if (Iteration::outerSteps % Iteration::nVisualSave == 0)
		{
			ONEFLOW::AddCmdToList( kVisualizationTaskName     );
			ONEFLOW::AddCmdToList( kDumpPressureCoeffTaskName );
			ONEFLOW::AddCmdToList( kDumpHeatfluxCoeffTaskName );
		}
	}

	if (startStrategy == 2)
	{
		if (Iteration::innerSteps % Iteration::nFieldSave == 0)
		{
			ONEFLOW::AddCmdToList( kDumpRestartTaskName );
		}
	}
	else
	{
		if (Iteration::outerSteps % Iteration::nFieldSave == 0)
		{
			ONEFLOW::AddCmdToList( kDumpRestartTaskName );
		}
	}

}

void INsFinalPostprocess( StringField & data )
{
    ONEFLOW::AddCmdToList( kDumpRestartTaskName       );
    ONEFLOW::AddCmdToList( kDumpAerodynamicTaskName   );
    ONEFLOW::AddCmdToList( kDumpPressureCoeffTaskName );
    ONEFLOW::AddCmdToList( kDumpHeatfluxCoeffTaskName );
    ONEFLOW::AddCmdToList( kVisualizationTaskName     );
}

void INsInitSolver( StringField & data )
{
    ONEFLOW::CommInterfaceData();
}

void IDumpHeatFluxCoeff( StringField & data )
{
    ;
}

void INsCalcTurb(StringField & data)
{
	;
}

void INsCalcHeat(StringField & data)
{
	;
}

SolverRegData insReg;
SolverRegData * GetINsReg()
{
    insReg.solverType = INC_NS_SOLVER;
    insReg.func = & RegisterINsFunc;
    insReg.solverName = "ins";
    insReg.baseKind = SOLVER_BASED;
    insReg.dataFlag = WITH_DATA;
    return & insReg;
}

REGISTER_REG_DATA( GetINsReg );

EndNameSpace
