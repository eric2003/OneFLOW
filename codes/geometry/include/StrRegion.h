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
#include "HXPointer.h"

BeginNameSpace( ONEFLOW )

class MyRegion
{
public:
    MyRegion();
    ~MyRegion();
public:
    IntField ijkmin, ijkmax;
public:
    int GetDirection();
};

typedef HXPointer< MyRegion > MyRegions;

class MyRRegion
{
public:
    MyRRegion();
    ~MyRRegion();
public:
    IntField idiv, jdiv, kdiv;
    MyRegions subregions;
    MyRegions refregions;
    MyRegions bcregions;
    MyRegions regions_nobc;
public:
    void CalcDiv( MyRegions & regions );
    void GenerateRegions( MyRegions & regions );
    void CollectNoSetBoundary();
    bool InBoundary( MyRegion * region );
    bool InRegion( MyRegion * r1, MyRegion * r2 );
    void AddRegion( MyRegion * region );
    void AddRefRegion( MyRegion * region );
    void AddRefRegion( MyRegions & regions );
    void AddBcRegion( MyRegion * region );
    void AddBcRegion( MyRegions & regions );
public:
    void Test();
    void Run();
};

class MyRegionFactory
{
public:
    MyRegionFactory();
    ~MyRegionFactory();
public:
    int ni, nj, nk;
    MyRegions refregions;
    MyRegions ref_bcregions;
    MyRegions bcregions;
public:
    void CreateRegion();
    void Create( int imin, int imax, int jmin, int jmax, int kmin, int kmax );
    void AddRefBcRegion( IntField & ijkMin, IntField & ijkMax );
    void AddBcRegion( MyRegions & bcregions_notset );
public:
    void Run();
    void CollectBcRegion( MyRegion * r, MyRegions & bcregions_collect );
};

EndNameSpace