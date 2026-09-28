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
#include <stdexcept>
#include <string>
#include <utility>
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
class ZgridMediator
{
public:
    ZgridMediator() = default;
    ~ZgridMediator() = default;

    ZgridMediator( const ZgridMediator & ) = delete;
    ZgridMediator & operator=( const ZgridMediator & ) = delete;
    ZgridMediator( ZgridMediator && ) noexcept = default;
    ZgridMediator & operator=( ZgridMediator && ) noexcept = default;

public:
    void AddGridMediator( GridMediator * gridMediator );
    void AddGridMediator( std::unique_ptr< GridMediator > gridMediator );

    [[nodiscard]] GridMediator * GetGridMediator( int iGridMediator ) const;
    [[nodiscard]] int GetSize() const;
    [[nodiscard]] std::string GetTargetFile() const;

public:
    void CreateSimple( int nZone );
    void ReadGrid();

private:
    std::vector< std::unique_ptr< GridMediator > > gm;
};

// Process-wide "current" GridMediator for legacy call paths that cannot
// take an explicit pointer. Prefer ScopedCurrentGridMediator in new code.
class GlobalGrid
{
public:
    GlobalGrid() = default;
    ~GlobalGrid() = default;

    // Non-owning. Caller must ensure lifetime exceeds all GetGrid uses.
    static void SetCurrentGridMediator( GridMediator * gridMediatorIn );

    [[nodiscard]] static GridMediator * GetCurrentGridMediator() noexcept;

    // Throws std::logic_error if no mediator is installed.
    [[nodiscard]] static Grid * GetGrid( int zoneId );

    // Historical public data member - prefer GetCurrentGridMediator().
    // Kept so existing TU that read GlobalGrid::gridMediator still link.
    static GridMediator * gridMediator;
};

// RAII: installs a current GridMediator for the enclosing scope and
// restores the previous one on destruction (including stack unwind).
class ScopedCurrentGridMediator
{
public:
    explicit ScopedCurrentGridMediator( GridMediator * next )
        : previous_( GlobalGrid::GetCurrentGridMediator() )
    {
        GlobalGrid::SetCurrentGridMediator( next );
    }

    ~ScopedCurrentGridMediator()
    {
        GlobalGrid::SetCurrentGridMediator( previous_ );
    }

    ScopedCurrentGridMediator( const ScopedCurrentGridMediator & ) = delete;
    ScopedCurrentGridMediator & operator=( const ScopedCurrentGridMediator & ) = delete;

private:
    GridMediator * previous_;
};

EndNameSpace
