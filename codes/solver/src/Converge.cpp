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

#include "Converge.h"
#include "Zone.h"
#include "ZoneState.h"
#include "DataBook.h"
#include "ActionState.h"
#include "SolverState.h"
#include "BasicParallel.h"

BeginNameSpace( ONEFLOW )

ConvergeTask::ConvergeTask()
{
    ;
}

ConvergeTask::~ConvergeTask()
{
    ;
}

void ConvergeTask::Run()
{
    ActionState::dataBook = this->dataBook;

    boolField.resize( ZoneState::nLocal, false );

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        if (  ! ZoneState::IsValidZone( zId ) ) continue;

        ZoneState::zid = zId;
        Solver * solver = SolverState::GetSolver();

    }

    this->CalcBool();
}

void ConvergeTask::CalcBool()
{
    this->flag = true;
    int nSize = boolField.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        this->flag &= boolField[ i ];
    }
    int s = this->flag;
    int t = -1;
    HXReduceInt( & s, & t, 1, PL_MIN );
    this->flag = t;
}

EndNameSpace