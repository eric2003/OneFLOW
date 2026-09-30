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
#include <string>

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
struct GridConfig;

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
    void GenerateGrid( const std::string & caseDir );
    void GenerateGrid( const GridConfig & config, const std::string & caseDir );
    void ReadCgnsGrid( const std::string & caseDir );
    void ReadCgnsGrid( const GridConfig & config, const std::string & caseDir );
    void DumpCgnsGrid( ZgridMediator & zgridMediator );
    void DumpUnsCgnsGrid( const std::string & caseDir );
    void DumpUnsCgnsGrid( const GridConfig & config, const std::string & caseDir );
public:
    void CommonToOneFlowGrid();
    void CommonToOneFlowGrid( const GridConfig & config );
    void CommonToStrGrid();
    void CommonToUnsGridTEST();
    void CommonToUnsGridTEST( const GridConfig & config );
    void ReadGridAndConvertToUnsCgnsZone();
    void ReadGridAndConvertToUnsCgnsZone( const GridConfig & config );
    void ProcessCgnsBases();
public:
    void CreateCgnsZone( ZgridMediator & zgridMediator );
    void PrepareCgnsZone( ZgridMediator & zgridMediator );
    CgnsZone * CreateSu2CgnsZone( Su2Grid & su2Grid );
    void Su2ToOneFlowGrid( Su2Grid & su2Grid );
public:
    void CgnsToOneFlowGrid();
    void CgnsToOneFlowGrid( const GridConfig & config );
    void ConvertStrCgns2UnsCgnsGrid();
};

void AddOneFlowGrid( Grids & grids, Grid * grid );
void GenerateLocalOneFlowGridFromSu2Grid( Su2Grid & su2Grid, Grids & grids );

#endif

EndNameSpace
