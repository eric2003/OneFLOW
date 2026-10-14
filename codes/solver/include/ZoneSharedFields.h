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
#include "HXArray.h"
#include "IFieldLookup.h"

BeginNameSpace( ONEFLOW )

class Grid;

// Zone-level shared physics fields (not owned by any single solver).
// Names are zone contracts, not SolverInfo slots.
struct ZoneSharedFieldNames
{
    static constexpr const char * ViscosityLaminar   = "visl";
    static constexpr const char * ViscosityTurbulent = "vist";
    static constexpr const char * TimeStep           = "timestep";
    static constexpr const char * Gamma              = "gama";
};

struct ZoneSharedFieldView
{
    MRField * visl     = nullptr;
    MRField * vist     = nullptr;
    MRField * timestep = nullptr;
    MRField * gama     = nullptr;

    // Injectable Bind -- unit-testable without DataBase.
    void Bind( Grid * grid, IFieldLookup & lookup );

    // Production: uses DataBaseFieldLookup (defined with DataBase adapter).
    void Bind( Grid * grid );
};

EndNameSpace
