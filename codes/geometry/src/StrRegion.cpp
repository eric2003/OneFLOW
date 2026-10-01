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

#include "StrRegion.h"
#include "HXStd.h"
#include "HXMath.h"
#include <iostream>
#include <memory>
#include <utility>


BeginNameSpace( ONEFLOW )

MyRegion::MyRegion()
{
}

MyRegion::~MyRegion()
{
}

int MyRegion::GetDirection()
{
    int dir = -1;
    for ( int i = 0; i < static_cast< int >( ijkmin.size() ); ++ i )
    {
        int var1 = ijkmin[ i ];
        int var2 = ijkmax[ i ];
        if ( var1 == var2 )
        {
            if ( var1 == 1 )
            {
                return 2 * i;
            }
            else
            {
                return 2 * i + 1;
            }
        }
    }
    return dir;
}

MyRRegion::MyRRegion()
{
}

void MyRRegion::AddRefRegion( MyRegion * region )
{
    this->refregions.push_back( region );
}

void MyRRegion::AddRefRegion( MyRegionViews & regions )
{
    for ( std::size_t i = 0; i < regions.size(); ++ i )
    {
        this->AddRefRegion( regions[ i ] );
    }
}

void MyRRegion::AddBcRegion( MyRegion * region )
{
    this->bcregions.push_back( region );
}

void MyRRegion::AddBcRegion( MyRegionViews & regions )
{
    for ( std::size_t i = 0; i < regions.size(); ++ i )
    {
        this->AddBcRegion( regions[ i ] );
    }
}

void MyRRegion::Test()
{
    MyRegion r1;
    r1.ijkmin.push_back( 1 );
    r1.ijkmin.push_back( 1 );
    r1.ijkmin.push_back( 1 );

    r1.ijkmax.push_back( 1 );
    r1.ijkmax.push_back( 5 );
    r1.ijkmax.push_back( 10 );

    MyRegion r2;
    r2.ijkmin.push_back( 1 );
    r2.ijkmin.push_back( 1 );
    r2.ijkmin.push_back( 1 );

    r2.ijkmax.push_back( 1 );
    r2.ijkmax.push_back( 3 );
    r2.ijkmax.push_back( 10 );

    this->AddRefRegion( & r1 );
    this->AddRefRegion( & r2 );

    this->AddBcRegion( & r2 );

    this->CalcDiv( refregions );

    this->GenerateRegions( subregions );

    this->CollectNoSetBoundary();
}

void MyRRegion::Run()
{
    this->CalcDiv( this->refregions );

    this->GenerateRegions( this->subregions );

    this->CollectNoSetBoundary();
}

void MyRRegion::CalcDiv( MyRegionViews & regions )
{
    IntSet idiv_set, jdiv_set, kdiv_set;
    const int nr = static_cast< int >( regions.size() );

    for ( int ir = 0; ir < nr; ++ ir )
    {
        MyRegion * region = regions[ ir ];
        idiv_set.insert( region->ijkmin[ 0 ] );
        idiv_set.insert( region->ijkmax[ 0 ] );

        jdiv_set.insert( region->ijkmin[ 1 ] );
        jdiv_set.insert( region->ijkmax[ 1 ] );

        kdiv_set.insert( region->ijkmin[ 2 ] );
        kdiv_set.insert( region->ijkmax[ 2 ] );
    }

    ONEFLOW::Set2Array( idiv_set, idiv );
    ONEFLOW::Set2Array( jdiv_set, jdiv );
    ONEFLOW::Set2Array( kdiv_set, kdiv );
}

void MyRRegion::GenerateRegions( MyRegions & regions )
{
    int ni = static_cast< int >( idiv.size() );
    int nj = static_cast< int >( jdiv.size() );
    int nk = static_cast< int >( kdiv.size() );

    int imin = 1;
    int imax = 1;

    int jmin = 1;
    int jmax = 1;

    int kmin = 1;
    int kmax = 1;

    int nni = MAX( ni - 1, 1 );
    int nnj = MAX( nj - 1, 1 );
    int nnk = MAX( nk - 1, 1 );

    int di = 1;
    int dj = 1;
    int dk = 1;

    if ( ni == 1 ) di = 0;
    if ( nj == 1 ) dj = 0;
    if ( nk == 1 ) dk = 0;

    for ( int k = 0; k < nnk; ++ k )
    {
        kmin = kdiv[ k ];
        kmax = kdiv[ k + dk ];

        for ( int j = 0; j < nnj; ++ j )
        {
            jmin = jdiv[ j ];
            jmax = jdiv[ j + dj ];
            for ( int i = 0; i < nni; ++ i )
            {
                imin = idiv[ i ];
                imax = idiv[ i + di ];

                auto region = std::make_unique< MyRegion >();
                region->ijkmin.push_back( imin );
                region->ijkmin.push_back( jmin );
                region->ijkmin.push_back( kmin );

                region->ijkmax.push_back( imax );
                region->ijkmax.push_back( jmax );
                region->ijkmax.push_back( kmax );
                regions.push_back( std::move( region ) );
            }
        }
    }
}

void MyRRegion::CollectNoSetBoundary()
{
    for ( std::size_t i = 0; i < subregions.size(); ++ i )
    {
        MyRegion * region = RegionAt( subregions, i );
        if ( ! this->InBoundary( region ) )
        {
            this->AddRegion( region );
        }
    }
}

void MyRRegion::AddRegion( MyRegion * region )
{
    this->regions_nobc.push_back( region );
}

bool MyRRegion::InBoundary( MyRegion * region )
{
    for ( std::size_t i = 0; i < bcregions.size(); ++ i )
    {
        MyRegion * bc_region = bcregions[ i ];
        if ( this->InRegion( region, bc_region ) )
        {
            return true;
        }
    }
    return false;
}

bool MyRRegion::InRegion( MyRegion * r1, MyRegion * r2 )
{
    if ( r1->ijkmin[ 0 ] < r2->ijkmin[ 0 ] ) return false;
    if ( r1->ijkmax[ 0 ] > r2->ijkmax[ 0 ] ) return false;

    if ( r1->ijkmin[ 1 ] < r2->ijkmin[ 1 ] ) return false;
    if ( r1->ijkmax[ 1 ] > r2->ijkmax[ 1 ] ) return false;

    if ( r1->ijkmin[ 2 ] < r2->ijkmin[ 2 ] ) return false;
    if ( r1->ijkmax[ 2 ] > r2->ijkmax[ 2 ] ) return false;
    return true;
}

MyRegionFactory::MyRegionFactory()
{
}

void MyRegionFactory::CreateRegion()
{
    this->Create( 1 , 1 , 1 , nj, 1 , nk );
    this->Create( ni, ni, 1 , nj, 1 , nk );

    this->Create( 1 , ni, 1 , 1 , 1 , nk );
    this->Create( 1 , ni, nj, nj, 1 , nk );

    this->Create( 1 , ni, 1 , nj, 1 , 1  );
    this->Create( 1 , ni, 1 , nj, nk, nk );
}

void MyRegionFactory::Create( int imin, int imax, int jmin, int jmax, int kmin, int kmax )
{
    auto r = std::make_unique< MyRegion >();
    r->ijkmin.push_back( imin );
    r->ijkmin.push_back( jmin );
    r->ijkmin.push_back( kmin );

    r->ijkmax.push_back( imax );
    r->ijkmax.push_back( jmax );
    r->ijkmax.push_back( kmax );

    this->refregions.push_back( std::move( r ) );
}

void MyRegionFactory::AddRefBcRegion( IntField & ijkMin, IntField & ijkMax )
{
    auto r = std::make_unique< MyRegion >();
    r->ijkmin = ijkMin;
    r->ijkmax = ijkMax;
    this->ref_bcregions.push_back( std::move( r ) );
}

void MyRegionFactory::AddBcRegion( MyRegionViews & bcregions_notset )
{
    for ( std::size_t i = 0; i < bcregions_notset.size(); ++ i )
    {
        auto r = std::make_unique< MyRegion >();
        MyRegion * rr = bcregions_notset[ i ];
        r->ijkmin = rr->ijkmin;
        r->ijkmax = rr->ijkmax;
        this->bcregions.push_back( std::move( r ) );
    }
}

void MyRegionFactory::Run()
{
    const int nFaces = static_cast< int >( this->refregions.size() );
    for ( int i = 0; i < nFaces; ++ i )
    {
        MyRegion * r = RegionAt( this->refregions, i );
        MyRegionViews bcregions_collect;
        this->CollectBcRegion( r, bcregions_collect );
        MyRRegion myrr;
        myrr.AddRefRegion( r );
        myrr.AddRefRegion( bcregions_collect );
        myrr.AddBcRegion( bcregions_collect );
        myrr.Run();
        AddBcRegion( myrr.regions_nobc );
    }
}

void MyRegionFactory::CollectBcRegion( MyRegion * r, MyRegionViews & bcregions_collect )
{
    for ( std::size_t i = 0; i < this->ref_bcregions.size(); ++ i )
    {
        MyRegion * rbc = RegionAt( this->ref_bcregions, i );
        if ( rbc->GetDirection() == r->GetDirection() )
        {
            bcregions_collect.push_back( rbc );
        }
    }
}

EndNameSpace
