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
#include "GridHandles.h"
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

BeginNameSpace( ONEFLOW )

class Grid;
struct GridConfig;

// Holds one zone-group of grids plus the paths / format used to load them.
// Historical name "Mediator" is kept for source compatibility; think of it as
// a GridBundle / GridDocument.
class GridMediator
{
public:
    GridMediator() = default;
    ~GridMediator() = default;

public:
    Grids gridVector;                 // owned grids (one entry per zone)
    int numberOfZones{ 0 };           // usually equals gridVector.size()
    int readGridType{ 0 };            // legacy flag; prefer gridType string / GridFileType
    std::string gridFile;             // source mesh path
    std::string bcFile;               // boundary condition path
    std::string targetFile;           // conversion output path
    std::string gridType;             // format token: plot3d, gridgen, ...
    std::string caseDir;              // explicit case root for grid file IO

public:
    void ReadGrid();
    void ReadGridgen();
    void ReadPlot3D();
    void ReadPlot3DCoor();
    void AddDefaultName();
};

// Owns a list of GridMediator instances (RAII).
// Name kept as ZgridMediator for compatibility; preferred mental model:
// "ZoneGridMediators" - a container of per-zone-group mediators.
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
    // --- modern container-style API (prefer these in new code) ---
    void add( std::unique_ptr< GridMediator > mediator );
    void add( GridMediator * mediator ); // takes ownership of a raw new'd pointer

    [[nodiscard]] GridMediator * at( int index );
    [[nodiscard]] const GridMediator * at( int index ) const;
    [[nodiscard]] int size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::string targetFile() const;

    // --- historical names (thin wrappers; keep call sites compiling) ---
    void AddGridMediator( GridMediator * gridMediator ) { add( gridMediator ); }
    void AddGridMediator( std::unique_ptr< GridMediator > gridMediator )
    {
        add( std::move( gridMediator ) );
    }

    [[nodiscard]] GridMediator * GetGridMediator( int iGridMediator )
    {
        return at( iGridMediator );
    }

    [[nodiscard]] const GridMediator * GetGridMediator( int iGridMediator ) const
    {
        return at( iGridMediator );
    }

    [[nodiscard]] int GetSize() const noexcept { return size(); }

    [[nodiscard]] std::string GetTargetFile() const { return targetFile(); }

public:
    void CreateSimple( int nZone );
    void ReadGrid();
    void ReadGrid( const GridConfig & config );

private:
    std::vector< std::unique_ptr< GridMediator > > mediators_;
};

EndNameSpace
