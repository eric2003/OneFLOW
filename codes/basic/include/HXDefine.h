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
#include "HXVector.h"
#include "HXType.h"
#include <set>
#include <string>

BeginNameSpace( ONEFLOW )

using RealField = HXVector< Real >;
using RealField2D = HXVector< RealField >;
using RealField3D = HXVector< RealField2D >;
using IntField = HXVector< int >;
using LinkField = HXVector< IntField >;
using StringField = HXVector< std::string >;
using BoolField = HXVector< bool >;

using IntSet = std::set< int >;

using VoidFunc = void( * )();

EndNameSpace
