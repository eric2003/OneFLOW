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
#include <cstddef>
#include <memory>
#include <vector>

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

// Owning region list.
using MyRegions = std::vector< std::unique_ptr< MyRegion > >;
// Non-owning observers / aliases (must not outlive the owners).
using MyRegionViews = std::vector< MyRegion * >;

[[nodiscard]] inline MyRegion * RegionAt( MyRegions & regions, std::size_t i )
{
    return regions[ i ].get();
}

[[nodiscard]] inline MyRegion * RegionAt( const MyRegions & regions, std::size_t i )
{
    return regions[ i ].get();
}

class MyRRegion
{
public:
    MyRRegion();
    ~MyRRegion() = default;
public:
    IntField idiv, jdiv, kdiv;
    MyRegions subregions;          // owned subdivisions
    MyRegionViews refregions;      // non-owning reference regions
    MyRegionViews bcregions;       // non-owning BC regions
    MyRegionViews regions_nobc;    // non-owning aliases into subregions
public:
    void CalcDiv( MyRegionViews & regions );
    void GenerateRegions( MyRegions & regions );
    void CollectNoSetBoundary();
    bool InBoundary( MyRegion * region );
    bool InRegion( MyRegion * r1, MyRegion * r2 );
    void AddRegion( MyRegion * region );
    void AddRefRegion( MyRegion * region );
    void AddRefRegion( MyRegionViews & regions );
    void AddBcRegion( MyRegion * region );
    void AddBcRegion( MyRegionViews & regions );
public:
    void Test();
    void Run();
};

class MyRegionFactory
{
public:
    MyRegionFactory();
    ~MyRegionFactory() = default;
public:
    int ni, nj, nk;
    MyRegions refregions;       // owned
    MyRegions ref_bcregions;    // owned
    MyRegions bcregions;        // owned
public:
    void CreateRegion();
    void Create( int imin, int imax, int jmin, int jmax, int kmin, int kmax );
    void AddRefBcRegion( IntField & ijkMin, IntField & ijkMax );
    void AddBcRegion( MyRegionViews & bcregions_notset );
public:
    void Run();
    void CollectBcRegion( MyRegion * r, MyRegionViews & bcregions_collect );
};

EndNameSpace
