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
#include "GridDef.h"
#include <memory>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

class Grid;

class GridMediator
{
public:
    GridMediator() = default;
    ~GridMediator() = default;

public:
    Grids gridVector;
    int numberOfZones{ 0 };
    int readGridType{ 0 };
    std::string gridFile;
    std::string bcFile;
    std::string targetFile;
    std::string gridType;

public:
    void ReadGrid();
    void ReadGridgen();
    void ReadPlot3D();
    void ReadPlot3DCoor();
    void AddDefaultName();
};

// Owns a collection of GridMediator instances (RAII).
// Historical SetDeleteFlag is kept as a no-op for source compatibility;
// ownership is always active for mediators added via AddGridMediator /
// CreateSimple / ReadGrid.
class ZgridMediator
{
public:
    ZgridMediator() = default;
    ~ZgridMediator() = default;

    // Non-copyable (owns unique resources); movable.
    ZgridMediator( const ZgridMediator & ) = delete;
    ZgridMediator & operator=( const ZgridMediator & ) = delete;
    ZgridMediator( ZgridMediator && ) noexcept = default;
    ZgridMediator & operator=( ZgridMediator && ) noexcept = default;

public:
    // Takes ownership of a heap-allocated GridMediator.
    // Prefer the unique_ptr overload in new code.
    void AddGridMediator( GridMediator * gridMediator );
    void AddGridMediator( std::unique_ptr< GridMediator > gridMediator );

    // Non-owning observer; valid while this ZgridMediator lives.
    [[nodiscard]] GridMediator * GetGridMediator( int iGridMediator ) const;
    [[nodiscard]] int GetSize() const;
    [[nodiscard]] std::string GetTargetFile() const;

    // Historical API: no longer needed. Ownership is always enabled.
    // Kept so existing call sites compile without change.
    void SetDeleteFlag( bool /*flag*/ ) {}

public:
    void CreateSimple( int nZone );
    void ReadGrid();

private:
    std::vector< std::unique_ptr< GridMediator > > gm;
};

class GlobalGrid
{
public:
    GlobalGrid() = default;
    ~GlobalGrid() = default;

public:
    static GridMediator * gridMediator;
    static Grid * GetGrid( int zoneId );
    static void SetCurrentGridMediator( GridMediator * gridMediatorIn );
};

EndNameSpace
