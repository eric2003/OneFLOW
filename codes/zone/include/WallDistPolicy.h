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


#pragma once
#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

enum class FlowWallDistAction
{
    SkipAfterAlloc,  // inviscid / laminar NS: allocate only
    Load,            // read existing wall distance
    Create           // compute and write wall distance
};

// Pure decision ¡ª no I/O, no tasks. Safe for unit tests.
inline FlowWallDistAction DecideFlowWallDistAction(
    int vismodel,
    int startStrategy,
    int ireadwdst )
{
    if ( vismodel <= 1 )
    {
        return FlowWallDistAction::SkipAfterAlloc;
    }
    if ( startStrategy > 0 )
    {
        return FlowWallDistAction::Load;
    }
    if ( ireadwdst == 0 )
    {
        return FlowWallDistAction::Create;
    }
    return FlowWallDistAction::Load;
}

EndNameSpace
