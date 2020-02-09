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
#include "SolverTaskReg.h"
#include "SolverRegData.h"
#include "Register.h"
#include "SolverDef.h"
#include "Category.h"
#include "RegisterUtil.h"
#include "SolverInfo.h"
#include "SolverName.h"

BeginNameSpace( ONEFLOW )

void RegisterSolverVarMap( int sTid )
{
    VarNameFactory::AddVarNameSolver( sTid, ONEFLOW::INTERFACE_DATA );
    VarNameFactory::AddVarNameSolver( sTid, ONEFLOW::INTERFACE_DQ_DATA );
    VarNameFactory::AddVarNameSolver( sTid, ONEFLOW::INTERFACE_GRADIENT_DATA );
    VarNameFactory::AddVarNameSolver( sTid, ONEFLOW::INTERFACE_OVERSET_DATA );
    SolverInfoFactory::AddSolverInfo( sTid );
}

void RegisterSolverTask( SolverRegData * regData )
{
    int sTid = regData->sTid;
    string &solverName = regData->solverName;
    VoidFunc func = regData->func;
    int baseKind = regData->baseKind;
    int dataFlag = regData->dataFlag;

    RegisterFactory::AddMRegister( sTid );
    Category::AddCategory( sTid, baseKind );

    if ( dataFlag == ONEFLOW::WITH_DATA )
    {
        RegisterSolverVarMap( sTid );
    }

    RegisterSolverFunc( sTid, solverName, func );
}

void RegisterSolverFunc( int sTid, const string & solverName, VoidFunc func )
{
    func();
    MRegister * mRegister = RegisterFactory::GetMRegister( sTid );
    SetSolverFileNames( mRegister, solverName );
    mRegister->RegisterAll();
}

void FreeSolverTask()
{
    RegisterFactory::FreeMRegister();
    Category::Free();
    VarNameFactory::FreeVarNameSolver();
    SolverInfoFactory::Free();
}

void SetSolverFileNames( MRegister * mRegister, const string & solverName )
{
    StringField fileNameList;
    ONEFLOW::GetSolverFileNames( solverName, fileNameList );
    mRegister->SetSolverFileNames( fileNameList );
}


EndNameSpace