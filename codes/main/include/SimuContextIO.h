/*---------------------------------------------------------------------------*\\
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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#pragma once

#include "SimuContextParameter.h"

BeginNameSpace( ONEFLOW )

namespace SimuContextDetail
{

// Apply a validated case-I/O snapshot to the legacy project-path state and
// write the diagnostic files. The caller must bind the case root first; the
// production setup calls this after ReadControlInfo(caseDir).
void ApplyCaseIOSettings( const CaseIOSettings & settings );

}

EndNameSpace
