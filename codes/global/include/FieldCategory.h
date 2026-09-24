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

// Defines the grid types to which a field allocation is applicable.
// All means no grid-type restriction; it is not a concrete grid type.
enum class FieldApplicability
{
    // Field is applicable to all supported grid types.
    All,

    // Field is applicable only to structured grids.
    Structured,

    // Field is applicable only to unstructured grids.
    Unstructured
};

enum class FieldLocation
{
    // Field is stored in the inner region.
    Inner,

    // Field is stored on faces.
    Face,

    // Field is stored in the boundary region.
    Boundary
};

EndNameSpace