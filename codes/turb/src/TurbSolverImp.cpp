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

#include "TurbSolverImp.h"
#include "SolverImp.h"
#include "SolverDef.h"
#include "TurbRhs.h"
#include "CmxTask.h"
#include "Iteration.h"
#include "SolverRegister.h"

BeginNameSpace( ONEFLOW )

void SetTurbFunc()
{
    REGISTER_DATA_CLASS( TurbInitFinal );
    REGISTER_DATA_CLASS( TurbUpdateResiduals );
    REGISTER_DATA_CLASS( TurbImplicitMethod );
    REGISTER_DATA_CLASS( TurbPostprocess );
    REGISTER_DATA_CLASS( TurbFinalPostprocess );
    REGISTER_DATA_CLASS( TurbInitSolver );
    REGISTER_DATA_CLASS( TurbCalcBoundary );
}

void TurbInitFinal( StringField & data )
{
    CalcTurbulentViscosity();

    TurbCalcBc();
}

void TurbCalcBoundary( StringField & data )
{
    CalcTurbulentViscosity();

    TurbCalcBc();
}

void TurbUpdateResiduals( StringField & data )
{
    TurbRhs * rhs = new TurbRhs();
    rhs->CalcRHS();
    delete rhs;
}

void TurbImplicitMethod( StringField & data )
{
    ;
}

void TurbPostprocess( StringField & data )
{
    CommInterfaceData();

    //残差的求解、输出，需要经过逻辑判断
    if ( Iteration::ResOk() )
    {
        ONEFLOW::AddCmdToList( "DUMP_RESIDUAL" );
    }

    if ( ! Iteration::InnerOk() ) return;

    ONEFLOW::AddCmdToList( "UPDATE_UNSTEADY_FLOW" );

    if ( Iteration::outerSteps % Iteration::nFieldSave == 0 )
    {
        ONEFLOW::AddCmdToList( "DUMP_RESTART" );
    }
}

void TurbFinalPostprocess( StringField & data )
{
    ONEFLOW::AddCmdToList( "DUMP_RESTART" );
}

void TurbInitSolver( StringField & data )
{
    ONEFLOW::CommInterfaceData();
}


SolverRegData turbReg;
SolverRegData * GetTurbReg()
{
    turbReg.sTid = TURB_SOLVER;
    turbReg.func = & SetTurbFunc;
    turbReg.solverName = "turb";
    turbReg.baseKind = SOLVER_BASED;
    turbReg.dataFlag = WITH_DATA;
    return & turbReg;
}

REGISTER_REG_DATA( GetTurbReg );

EndNameSpace