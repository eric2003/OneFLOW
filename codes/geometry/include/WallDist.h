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

#include "HXClone.h"
#include "Point.h"
#include "Task.h"
#include <memory>

BeginNameSpace( ONEFLOW )

DEFINE_DATA_CLASS( FillWallStructTask );
DEFINE_DATA_CLASS( FillWallStruct );
DEFINE_DATA_CLASS( CalcWallDist );

class WallStructure
{
public:
    using PointType  = Point< Real >;
    using PointField = HXVector< PointType >;
    using PointLink  = HXVector< PointField >;

    PointField fc;
    PointLink  fv;
};

// Task that aggregates wall geometry across zones into the shared
// WallStructure storage (see GetWallStructure()).
class CFillWallStructTaskImp : public Task
{
public:
    CFillWallStructTaskImp() = default;
    ~CFillWallStructTaskImp() override = default;

    void Run() override;
    void Create();
    void FillWall();
};

// Shared wall-structure storage used between FILL_WALL_STRUCT and
// CALC_WALL_DIST tasks. Owned by unique_ptr; never a bare new/delete pair.
[[nodiscard]] WallStructure * GetWallStructure() noexcept;
void ResetWallStructure();           // release storage (replaces FreeWallStruct body)
void EnsureWallStructure();         // create if missing

void SetWallTask();
void FreeWallStruct();              // kept for existing call sites -> ResetWallStructure()

Real CalcPoint2FaceDist( WallStructure::PointType node,
                         WallStructure::PointField & fvList );

EndNameSpace
