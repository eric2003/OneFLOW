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

#include "RestartTaskReg.h"
#include "InterField.h"
#include "ActionState.h"
#include "HXMath.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Restart.h"
#include "Unsteady.h"
#include "UnsteadyImp.h"
#include "Update.h"
#include "FieldWrap.h"
#include "FieldAlloc.h"
#include "CmxTask.h"
#include "DataBase.h"
#include "DataBook.h"
#include "Lusgs.h"
#include "Lhs.h"
#include "FieldImp.h"
#include "FieldWrap.h"
#include "SolverState.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "InterFace.h"
#include "RegisterUtil.h"
#include "FieldRecord.h"
#include "UVisualize.h"
#include "UResidual.h"
#include "UNsCom.h"
#include "TaskRegister.h"
#include "UINsCom.h"
#include <map>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterRestartTask )

void RegisterRestartTask()
{
    REGISTER_DATA_CLASS( InitFirst );
    REGISTER_DATA_CLASS( ReadRestart );
    REGISTER_DATA_CLASS( DumpRestart );
    REGISTER_DATA_CLASS( InitRestart );
    REGISTER_DATA_CLASS( InitFlowField );

	REGISTER_DATA_CLASS( ReadinsRestart );
	REGISTER_DATA_CLASS( DumpinsRestart );
	REGISTER_DATA_CLASS( InitinsRestart );
}

void InitFirst( StringField & data )
{
    int sTid = SolverState::tid;

    const string & basicString = data[ 0 ];

    FieldAlloc::AllocateAllFields( sTid, basicString );
}

void ReadRestart( StringField & data )
{
    int sTid = SolverState::tid;

    Restart * restart = CreateRestart( sTid );
    restart->Read( sTid );
    delete restart;
}

void DumpRestart( StringField & data )
{
    int sTid = SolverState::tid;

    Restart * restart = CreateRestart( sTid );
    restart->Dump( sTid );
    delete restart;
}

void DumpinsRestart(StringField & data)
{
	int sTid = SolverState::tid;

	Restart * restart = CreateRestart(sTid);
	restart->Dump(sTid);
	delete restart;
}

void InitRestart( StringField & data )
{
    int sTid = SolverState::tid;

    Restart * restart = CreateRestart( sTid );
    restart->InitRestart( sTid );
    delete restart;
}

void InitinsRestart( StringField & data )
{
	int sTid = SolverState::tid;

	Restart * restart = CreateRestart( sTid );
	restart->InitinsRestart( sTid );
	delete restart;
}

void ReadinsRestart(StringField & data)
{
	int sTid = SolverState::tid;

	Restart * restart = CreateRestart(sTid);
	restart->Read( sTid );
	delete restart;
}


void InitFlowField( StringField & data )
{
    ONEFLOW::AddCmdToList( "INIT_FIRST" );

    int startStrategy = ONEFLOW::GetDataValue< int >( "startStrategy" );

    if ( startStrategy == 0 )
    {
        ONEFLOW::AddCmdToList( "INIT_RESTART" );
    }

	else if (startStrategy == 1)
	{
		ONEFLOW::AddCmdToList("READ_RESTART");
	}

	else if (startStrategy == 2)
	{
		ONEFLOW::AddCmdToList("INIT_INSRESTART");
	}

	else if ( startStrategy == 3 )
    {
        ONEFLOW::AddCmdToList( "READ_INSRESTART" );
    }

    ONEFLOW::AddCmdToList( "INIT_FINAL" );
}


EndNameSpace