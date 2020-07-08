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

#include "NsSolverImp.h"
#include "SolverImp.h"
#include "UTimeStep.h"
#include "Zone.h"
#include "GridState.h"
#include "UNsVisFlux.h"
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

void RegisterNsFunc()
{
    REGISTER_DATA_CLASS( NsInitFinal );
    REGISTER_DATA_CLASS( NsVisual );
    REGISTER_DATA_CLASS( NsCalcTimeStep );
    REGISTER_DATA_CLASS( NsUpdateResiduals );
    REGISTER_DATA_CLASS( NsImplicitMethod );
    REGISTER_DATA_CLASS( NsPostprocess );
    REGISTER_DATA_CLASS( NsFinalPostprocess );
    REGISTER_DATA_CLASS( NsInitSolver );
    REGISTER_DATA_CLASS( NsCalcBoundary );
    REGISTER_DATA_CLASS( DumpHeatFluxCoeff );
    SetPlateTask();
    SetTurbPlateTask();
}

void NsInitFinal( StringField & data )
{
    NsCalcGamaT( F_INNER );
    CalcLaminarViscosity( F_INNER );
    NsCalcBc();
    NsCalcGamaT( F_GHOST );
    CalcLaminarViscosity( F_GHOST );

    Grid * grid = Zone::GetGrid();

    if ( Zone::GetCGrid( grid ) )
    {
        //RestrictAllQ( NS_SOLVER, FLOW_FIELD_INDEX );

        GridState::gridLevel += 1;

        NsCalcBc();

        GridState::gridLevel -= 1;
    }
}

void NsVisual( StringField & data )
{
    ;
}

void NsCalcBoundary( StringField & data )
{
    NsCalcGamaT( F_INNER );
    CalcLaminarViscosity( F_INNER );
    NsCalcBc();
    NsCalcGamaT( F_GHOST );
    CalcLaminarViscosity( F_GHOST );
}

void NsCalcTimeStep( StringField & data )
{
    UTimeStep * uTimeStep = new UTimeStep();
    uTimeStep->CalcTimeStep();
    delete uTimeStep;
}

void NsUpdateResiduals( StringField & data )
{
    Rhs * rhs = new Rhs();
    rhs->UpdateNsResiduals();
    delete rhs;
}

void NsImplicitMethod( StringField & data )
{
    ;
}

void NsPostprocess( StringField & data )
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
        if ( ctrl.idump == 1 )
        {
            ONEFLOW::AddCmdToList( "DUMP_LAMINAR_PLATE" );
        }
        else if ( ctrl.idump == 2 )
        {
            ONEFLOW::AddCmdToList( "DUMP_TURB_PLATE" );
        }
    }

    if ( Iteration::outerSteps % Iteration::nFieldSave == 0 )
    {
        ONEFLOW::AddCmdToList( "DUMP_RESTART" );
    }

}

void DumpHeatFluxCoeff( StringField & data )
{
    ;
}

void NsFinalPostprocess( StringField & data )
{
    ONEFLOW::AddCmdToList( "DUMP_RESTART"        );
    ONEFLOW::AddCmdToList( "DUMP_AERODYNAMIC"    );
    ONEFLOW::AddCmdToList( "DUMP_PRESSURE_COEFF" );
    ONEFLOW::AddCmdToList( "DUMP_HEATFLUX_COEFF" );
    ONEFLOW::AddCmdToList( "VISUALIZATION"       );
}

void NsInitSolver( StringField & data )
{
    ONEFLOW::CommInterfaceData();
}

SolverRegData nsReg;
SolverRegData * GetNsReg()
{
    nsReg.sTid = NS_SOLVER;
    nsReg.func = & RegisterNsFunc;
    nsReg.solverName = "ns";
    nsReg.baseKind = SOLVER_BASED;
    nsReg.dataFlag = WITH_DATA;
    return & nsReg;
}

REGISTER_REG_DATA( GetNsReg );

EndNameSpace