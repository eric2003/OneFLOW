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

#include "Update.h"
#include "NsUpdate.h"
#include "INsUpdate.h"
#include "TurbUpdate.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Task.h"
#include "TaskState.h"
#include "FieldWrap.h"

BeginNameSpace( ONEFLOW )

Update::Update()
{
    q = 0;
    dq = 0;
}

Update::~Update()
{
    delete q;
    delete dq;
}

Update * CreateUpdate( int sTid )
{
    if ( sTid == NS_SOLVER )
    {
        return CreateNsUpdate();
    }
    else if ( sTid == INC_NS_SOLVER )
    {
        return CreateINsUpdate();
    }
    else if ( sTid == TURB_SOLVER )
    {
        return CreateTurbUpdate();
    }
    return 0;
}

void GetUpdateField( int sTid, FieldWrap *q, FieldWrap *dq )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );

    if ( TaskState::task->taskName == "UPDATE_FLOWFIELD_LUSGS" )
    {
        string & qFieldString  = solverInfo->implicitString[ 0 ];
        string & dQFieldString = solverInfo->implicitString[ 1 ];
        q  = FieldHome::GetFieldWrap( qFieldString  );
        dq = FieldHome::GetFieldWrap( dQFieldString );
    }
	else if (TaskState::task->taskName == "UPDATE_FLOWFIELD_SIMPLE")
	{
		string & qFieldString = solverInfo->implicitString[0];
		string & dQFieldString = solverInfo->implicitString[1];
		q = FieldHome::GetFieldWrap(qFieldString);
		dq = FieldHome::GetFieldWrap(dQFieldString);
	}
    else
    {
        q  = FieldHome::GetFieldWrap( FIELD_FLOW );
        dq = FieldHome::GetFieldWrap( solverInfo->residualName );
    }
}

EndNameSpace