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

#include "TurbSolver.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "TurbCtrl.h"
#include "TurbCom.h"
#include "DataBase.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_SOLVER( TurbSolver )

bool TurbSolver::initFlag = false;

TurbSolver::TurbSolver()
{
}

TurbSolver::~TurbSolver()
{
}

void TurbSolver::StaticInit()
{
    if ( TurbSolver::initFlag ) return;
    TurbSolver::initFlag = true;

    this->sTid = ONEFLOW::TURB_SOLVER;
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( this->sTid );
    solverInfo->nEqu  = 1;
    solverInfo->nTEqu = 1;

    if ( vis_model.vismodel == 3 )
    {
        solverInfo->nEqu  = 1;
        solverInfo->nTEqu = 1;
    }
    else if ( vis_model.vismodel == 4 )
    {
        solverInfo->nEqu  = 2;
        solverInfo->nTEqu = 2;
    }

    solverInfo->registerInterface = 0;
    solverInfo->residualName = "turbres";
    solverInfo->resFileName = GetDataValue< string >( "turbresFile" );
    solverInfo->gradString.push_back( "turbq"    );
    solverInfo->gradString.push_back( "turbdqdx" );
    solverInfo->gradString.push_back( "turbdqdy" );
    solverInfo->gradString.push_back( "turbdqdz" );

    solverInfo->implicitString.push_back( "turbq"  );
    solverInfo->implicitString.push_back( "turbdq" );

    int nTurbEqu = solverInfo->nEqu;
    SetDataInt( "nTurbEqu", nTurbEqu );
    turb_ctrl.Init();
    turbcom.Init();
}

EndNameSpace