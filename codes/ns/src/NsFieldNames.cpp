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

#include "NsFieldNames.h"

BeginNameSpace( ONEFLOW )

SolverLogisticsNames NsSolverFieldLogisticsNames()
{
    return MakeSolverLogisticsNames(
        "q", "dqdx", "dqdy", "dqdz", "dq", "res",
        "bcdqdx", "bcdqdy", "bcdqdz",
        "limiter", "bc_q" );
}

UnsteadyFieldNames NsSolverUnsteadyFieldNames()
{
    StringField stageQ;
    stageQ.push_back( "q1" );
    stageQ.push_back( "q2" );
    StringField stageR;
    stageR.push_back( "res1" );
    stageR.push_back( "res2" );
    return MakeUnsteadyFieldNames( "q", stageQ, "res", stageR );
}

EndNameSpace
