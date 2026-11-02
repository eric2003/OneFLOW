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
#include "System.h"
#include "DimensionImp.h"
#include "SolverRegister.h"
#include "SolverTaskReg.h"
#include "SolverDef.h"
#include "TaskRegister.h"
#include "MessageMapLoader.h"

BeginNameSpace( ONEFLOW )

void ConstructSystemMap()
{
    ONEFLOW::SetDimension();

    TaskRegister::Run();

    CreateSysMap();

    CreateMsgMap();

    // Rebuild solver registrations for each case. ConstructSystemMap() is
    // case-scoped in the multi-case execution path, so stale registration
    // data must not accumulate across cases.
    FreeSolverTask();
    SolverRegister::Run();
}

EndNameSpace
