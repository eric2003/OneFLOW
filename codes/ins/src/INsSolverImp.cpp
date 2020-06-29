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

#include "INsSolverImp.h"
#include "SolverImp.h"
#include "UITimeStep.h"
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
}

void INsInitFinal( StringField & data )
{
    INsCalcGamaT( F_INNER );
    ICalcLaminarViscosity( F_INNER );
    INsCalcBc();
    INsCalcGamaT( F_GHOST );
    ICalcLaminarViscosity( F_GHOST );
    Grid * grid = Zone::GetGrid();

    if ( Zone::GetCGrid( grid ) )
    {
        //RestrictAllQ( NS_SOLVER, FLOW_FIELD_INDEX );

        GridState::gridLevel += 1;

        NsCalcBc();

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
    ICalcLaminarViscosity( F_INNER );
    INsCalcBc();
    INsCalcGamaT( F_GHOST );
    ICalcLaminarViscosity( F_GHOST );
}

void INsCalcTimeStep( StringField & data )
{
    UITimeStep * uTimeStep = new UITimeStep();
    uTimeStep->CalcTimeStep();
    delete uTimeStep;
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

    if ( ! Iteration::InnerOk() ) return;

    ONEFLOW::AddCmdToList( "UPDATE_UNSTEADY_FLOW" );

    if ( Iteration::outerSteps % Iteration::nVisualSave == 0 )
    {
        ONEFLOW::AddCmdToList( "VISUALIZATION"       );
        ONEFLOW::AddCmdToList( "DUMP_PRESSURE_COEFF" );
        ONEFLOW::AddCmdToList( "DUMP_HEATFLUX_COEFF" );
        //if ( ctrl.idump == 1 )
        //{
        //    ONEFLOW::AddCmdToList( "DUMP_LAMINAR_PLATE" );
        //}
        //else if ( ctrl.idump == 2 )
        //{
        //    ONEFLOW::AddCmdToList( "DUMP_TURB_PLATE" );
        //}
    }

    if ( Iteration::outerSteps % Iteration::nFieldSave == 0 )
    {
        ONEFLOW::AddCmdToList( "DUMP_RESTART" );
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