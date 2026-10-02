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

#include "Update.h"
#include <memory>
#include "NsUpdate.h"
#include "INsUpdate.h"
#include "TurbUpdate.h"
#include "CmxTaskNames.h"
#include "SolverDef.h"
#include "SolverInfo.h"
#include "Task.h"
#include "TaskState.h"
#include "FieldWrap.h"

BeginNameSpace( ONEFLOW )

Update::Update()
{
}

Update::~Update()
{
}

Update * CreateUpdate( int solverType )
{
    if ( solverType == NS_SOLVER )
    {
        return CreateNsUpdate();
    }
    else if ( solverType == INC_NS_SOLVER )
    {
        return CreateINsUpdate();
    }
    else if ( solverType == TURB_SOLVER )
    {
        return CreateTurbUpdate();
    }
    return 0;
}

void GetUpdateField(
    int solverType,
    std::unique_ptr<FieldWrap> & q,
    std::unique_ptr<FieldWrap> & dq )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );

    if ( TaskState::task->taskName == kUpdateFlowFieldLusgsTaskName )
    {
        std::string & qFieldString  = solverInfo->implicitString[ 0 ];
        std::string & dQFieldString = solverInfo->implicitString[ 1 ];
        q.reset( FieldHome::GetFieldWrap( qFieldString  ) );
        dq.reset( FieldHome::GetFieldWrap( dQFieldString ) );
    }
    else
    {
        // FIELD_FLOW wrap is owned by BgField; only borrow its MRField via a
        // non-owning FieldWrap that Update uniquely owns.
        FieldWrap * bgFlow = FieldHome::GetFieldWrap( FIELD_FLOW );
        auto flowWrap = std::make_unique<FieldWrap>();
        flowWrap->SetUnsField( bgFlow->GetUnsField(), false );
        q = std::move( flowWrap );

        // residualName path builds a fresh non-owning wrapper; Update owns it.
        dq.reset( FieldHome::GetFieldWrap( solverInfo->residualName ) );
    }
}

EndNameSpace
