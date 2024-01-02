/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class PointFactory;
class CgnsZone;
class CgnsZbase;
class ElemFeature;
class FaceSolver;
class Grid;
class UnsGrid;
class CgnsSection;

int OneFlow2CgnsZoneType( int zoneType );
int Cgns2OneFlowZoneType( int zoneType );

class GridElem
{
public:
    GridElem( HXVector< CgnsZone * > & cgnsZones, int iZone );
    ~GridElem();
public:
    ElemFeature  * elem_feature;
    PointFactory * point_factory;
    FaceSolver   * face_solver;
    HXVector< CgnsZone * > cgnsZones;
    Grid * grid;
    bool delFlag;
    Real minLen, maxLen;
public:
    CgnsZone * GetCgnsZone( int iZone );
    int GetNZones();
public:
    void CreateGrid( HXVector< CgnsZone * > cgnsZones, int iZone );
    void PrepareUnsCalcGrid();
    void PrepareUnsCalcGridNormal();
    void InitCgnsElements();
    void ScanBcFace();
    void GenerateCalcElement();
    void GenerateCalcGrid();
    void GenerateCalcGrid( Grid * grid );
    void CalcBoundaryType( UnsGrid * grid );
    void ReorderLink( UnsGrid * grid );
public:
    void PrepareUnsCalcGridPolyhedron();
    void ScanPolygonFace();
    void SetPolyhedronElementType( CgnsSection * cgnsSection );
};

class ZgridElem
{
public:
    ZgridElem( CgnsZbase * cgnsZbase );
    ~ZgridElem();
public:
    HXVector< GridElem * > data;
    CgnsZbase * cgnsZbase;
    Grids grids;
public:
    GridElem * GetGridElem( int iGridElem );
    void AddGridElem( GridElem * gridElem );
    void AddGridElem( HXVector< CgnsZone * > cgnsZones, int iZone );
public:
    void GenerateLocalOneFlowGrid( Grids & grids );
    void AllocateGridElem();
    void PrepareUnsCalcGrid();
    void GenerateCalcGrid();
    void GetGrids( Grids & grids );
};

EndNameSpace
