/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "AdtTree.h"
#include "GridDef.h"

BeginNameSpace( ONEFLOW )

class Grid;

typedef HXAdtTree< int, Real > AdtTree;
typedef HXAdtNode< int, Real > AdtNode;

class PointSearch
{
public:
    PointSearch();
    ~PointSearch();
protected:
    int id;
    AdtTree * coorTree;
    Real tolerance;
    RealField xCoor, yCoor, zCoor;
public:
    void Initialize( RealField & pmin, RealField & pmax, Real toleranceIn );
    void Initialize( Grid * grid );
    void InitializeSpecial( Grid * grid, Real toleranceIn );
    void Initialize( Grids & grids );
public:
    int GetNPoint() { return static_cast<int> (xCoor.size()); }
    int FindPoint( Real xm, Real ym, Real zm );
    int AddPoint( Real xm, Real ym, Real zm );
    void GetPoint( int id, Real & xm, Real & ym, Real & zm );
    Real GetTol() { return tolerance; }
protected:
    int AddPoint( RealField & coordinate );
    int FindPoint( RealField & coordinate );
public:
    void GetFaceCoorList( IntField & nodeId, RealField &xList, RealField &yList, RealField &zList );
};

void CreateStandardADT( RealField & ptmin, RealField & ptmax, AdtTree *& adtTree, Real & tolerance );
void CreateStandardADT( Grid * grid, AdtTree *& adtTree, Real & tolerance );
void CreateStandardADT( Grids & grids, AdtTree *& adtTree, Real & tolerance );
void CreateStandardADTByTolerance( Grids & grids, AdtTree *& adtTree, Real & tolerance );

void ShiftMinMaxBox( RealField & pmin, RealField & pmax, Real tolerance );
void GetGridsMinMaxDistance( Grids & grids, Real & mindis, Real & maxdis );
Real CalcGridTolerance( Grids & grids );
void GetBoundingBoxOfMultiZoneGrids( Grids & grids, RealField & pmin, RealField & pmax );

EndNameSpace