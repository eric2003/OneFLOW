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
#include "HXDefine.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class Grid;

[[nodiscard]] std::unique_ptr< Grid > CreateGridUnique( int gridType );
[[nodiscard]] std::unique_ptr< Grid > CreateUnsGridUnique();
[[nodiscard]] std::unique_ptr< Grid > CreateStrGridUnique();

// Compatibility: returns a raw pointer the caller must own.
// Legacy raw-pointer factories (no remaining in-tree call sites).
// Prefer Create*Unique and store in Grids / unique_ptr.
[[deprecated( "Use CreateGridUnique" )]]
Grid * CreateGrid( int gridType );
[[deprecated( "Use CreateUnsGridUnique" )]]
Grid * CreateUnsGrid();
[[deprecated( "Use CreateStrGridUnique" )]]
Grid * CreateStrGrid();

EndNameSpace
