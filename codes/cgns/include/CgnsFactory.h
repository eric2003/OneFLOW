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
#include "GridDef.h"
#include "HXCgns.h"
#include <memory> // Required for std::unique_ptr

BeginNameSpace( ONEFLOW )

class StrGrid;
class NodeMesh;
class Grid;
class Su2Grid;
class CgnsZbase;
class CgnsZone;
class GridElem;
class ZgridElem;
class GridMediator;
class ZgridMediator;

#ifdef ENABLE_CGNS

class CgnsFactory
{
public:
    CgnsFactory();
    ~CgnsFactory();

    // Rule of 5: Disable copying to prevent double-free
    CgnsFactory(const CgnsFactory&) = delete;
    CgnsFactory& operator=(const CgnsFactory&) = delete;

    // FIX: Declare move semantics here, but DO NOT use = default.
    // The implementation must be in the .cpp file where types are complete.
    CgnsFactory(CgnsFactory&&) noexcept;
    CgnsFactory& operator=(CgnsFactory&&) noexcept;

public:
    std::unique_ptr<CgnsZbase> cgnsZbase;
    std::unique_ptr<ZgridElem> zgridElem;
public:
    void GenerateGrid();
    void ReadCgnsGrid();
    void DumpCgnsGrid( ZgridMediator * zgridMediator );
    void DumpUnsCgnsGrid();
public:
    void CommonToOneFlowGrid();
    void CommonToStrGrid();
    void CommonToUnsGridTEST();
    void ReadGridAndConvertToUnsCgnsZone();
    void ProcessCgnsBases();
public:
    void CreateCgnsZone( ZgridMediator * zgridMediator );
    void PrepareCgnsZone( ZgridMediator * zgridMediator );
    CgnsZone * CreateSu2CgnsZone( Su2Grid* su2Grid );
    void Su2ToOneFlowGrid( Su2Grid* su2Grid );
public:
    void CgnsToOneFlowGrid();
    void ConvertStrCgns2UnsCgnsGrid();
};

void AddOneFlowGrid( Grids & grids, Grid * grid );
void GenerateLocalOneFlowGridFromSu2Grid( Su2Grid* su2Grid, Grids & grids );

#endif

EndNameSpace
