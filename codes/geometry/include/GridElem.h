/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#pragma once
#include "HXDefine.h"
#include "GridHandles.h"
// These members are stored by value, so their types must be complete here.
#include "ElemFeature.h"
#include "PointManager.h"
#include "FaceSolver.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class MeshPointManager;
class CgnsZone;
class CgnsZbase;
class ElemFeature;
class FaceSolver;
class Grid;
class UnsGrid;
class CgnsSection;
struct GridConfig;

int OneFlow2CgnsZoneType( int zoneType );
int Cgns2OneFlowZoneType( int zoneType );

class GridElem
{
public:
    GridElem( const HXVector< CgnsZone * > & cgnsZones );
    ~GridElem();
public:
    ElemFeature elem_feature;
    MeshPointManager point_factory;
    FaceSolver face_solver;
    HXVector< CgnsZone * > cgnsZones;
    Real minLen, maxLen;
public:
    CgnsZone * GetCgnsZone( int iZone );
    const CgnsZone * GetCgnsZone( int iZone ) const;
    int GetNZones() const;
public:
    void PrepareUnsCalcGrid();
    void PrepareUnsCalcGridNormal();
    void InitCgnsElements();
    void ScanBcFace();
    void GenerateCalcElement();
    [[nodiscard]] std::unique_ptr< UnsGrid > GenerateCalcGrid( int gridId );
    void GenerateCalcGrid( UnsGrid & grid );
    void CalcBoundaryType( UnsGrid & grid );
    void ReorderLink( UnsGrid & grid );
public:
    void PrepareUnsCalcGridPolyhedron();
    void ScanPolygonFace();
    void SetPolyhedronElementType( CgnsSection & cgnsSection );
};

class ZgridElem
{
public:
    ZgridElem( CgnsZbase * cgnsZbase );
    ~ZgridElem();
private:
    CgnsZbase * cgnsZbase;
public:
    [[nodiscard]] CgnsZbase * GetCgnsZbase() const;
    void RebindCgnsZbase( CgnsZbase * cgnsZbase ) noexcept;

    [[nodiscard]] Grids GenerateLocalOneFlowGrids();
    [[nodiscard]] Grids GenerateLocalOneFlowGrids( const GridConfig & config );
private:
    [[nodiscard]] HXVector< std::unique_ptr< GridElem > > CreateGridElements( bool multiBlock ) const;
    void PrepareUnsCalcGrid( HXVector< std::unique_ptr< GridElem > > & data ) const;
};

EndNameSpace
