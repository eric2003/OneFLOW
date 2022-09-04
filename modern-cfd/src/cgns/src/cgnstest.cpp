/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

#include "cgnstest.h"
#include "cgnslib.h"
#include "tools.h"
#include <string>
#include <vector>
#include <iostream>

#ifdef PRJ_ENABLE_CGNS

int Cgns_t::cell_dim = 3;
int Cgns_t::phys_dim = 3;
int Cgns_t::file_id = -1;
int Cgns_t::base_id = -1;
int Cgns_t::zone_id = -1;
int Cgns_t::coor_id = -1;
int Cgns_t::bc_id   = -1;

std::vector< CgnsBase * > cgns_base_array;

void Cgns_t::WriteBase( const std::string & baseName )
{
    cg_base_write( Cgns_t::file_id, baseName.c_str(), Cgns_t::cell_dim, Cgns_t::phys_dim, & Cgns_t::base_id );
    std::cout << " CGNS Base index = " << Cgns_t::base_id << "\n";
}

void Cgns_t::SetISize( cgsize_t *isize, std::vector<int> &dim_array )
{
    int n = 0;
    // vertex size
    for ( int i = 0; i < dim_array.size(); ++ i )
    {
        isize[ n ++ ] = dim_array[ i ];
    }
    // cell size
    for ( int i = 0; i < dim_array.size(); ++ i )
    {
        isize[ n ++ ] = std::max(dim_array[ i ] - 1, 1);
    }
    // boundary vertex size (always zero for structured grids)
    for ( int i = 0; i < dim_array.size(); ++ i )
    {
        isize[ n ++ ] = 0;
    }
}

void Cgns_t::SetStructuredIpnts( cgsize_t *ipnts, int ilo, int jlo, int klo, int ihi, int jhi, int khi )
{
    int n = 0;
    ipnts[ n ++ ] = ilo;
    ipnts[ n ++ ] = jlo;
    ipnts[ n ++ ] = klo;

    ipnts[ n ++ ] = ihi;
    ipnts[ n ++ ] = jhi;
    ipnts[ n ++ ] = khi;
}

void Cgns_t::SetStructuredIpnts( cgsize_t *ipnts, int ilo, int jlo, int ihi, int jhi )
{
    int n = 0;
    ipnts[ n ++ ] = ilo;
    ipnts[ n ++ ] = jlo;

    ipnts[ n ++ ] = ihi;
    ipnts[ n ++ ] = jhi;
}

void Cgns_t::SetStructuredIpnts( cgsize_t *ipnts, int ilo )
{
    int n = 0;
    ipnts[ n ++ ] = ilo;
    ipnts[ n ++ ] = ilo;
}

void Cgns_t::SetDimArray( std::vector<int> &dim_array, int ni, int nj, int nk )
{
    dim_array.resize( 0 );
    dim_array.push_back( ni );
    dim_array.push_back( nj );
    dim_array.push_back( nk );
}

void Cgns_t::SetDimArray( std::vector<int> &dim_array, int ni, int nj )
{
    dim_array.resize( 0 );
    dim_array.push_back( ni );
    dim_array.push_back( nj );
}

void Cgns_t::SetDimArray( std::vector<int> &dim_array, int ni )
{
    dim_array.resize( 0 );
    dim_array.push_back( ni );
}

void cgns_write_base_test()
{
    if ( cg_open( "test.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim=1;
    Cgns_t::phys_dim=3;

    Cgns_t::WriteBase( "baseA");
    Cgns_t::WriteBase( "baseB");

    cg_close( Cgns_t::file_id );
}

int cgnstest()
{
    int index_file = -1;
    if ( cg_open( "test.cgns", CG_MODE_WRITE, &index_file ) )
    {
        cg_error_exit();
    }
    std::string basename = "Base";

    int icelldim=3;
    int iphysdim=3;
    int index_base = -1;
    cg_base_write(index_file,basename.c_str(), icelldim, iphysdim, &index_base);
    cg_close(index_file);
    return 0;
}

void cgns_dump_grid( float * xcoor, int ni, const std::string & filename )
{
    if ( cg_open( filename.c_str(), CG_MODE_WRITE, &Cgns_t::file_id) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 1;
    Cgns_t::phys_dim = 3;
    cg_base_write( Cgns_t::file_id, "Base", Cgns_t::cell_dim, Cgns_t::phys_dim, &Cgns_t::base_id );
    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    Cgns_t::SetDimArray( dim_array, ni );
    Cgns_t::SetISize( isize, dim_array );
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, "Zone 1", isize, CGNS_ENUMV(Structured), &Cgns_t::zone_id );
    cgsize_t ipnts[ 6 ];
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1 );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "left", CGNS_ENUMV( BCInflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, ni );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "right", CGNS_ENUMV( BCOutflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    cg_coord_write ( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, CGNS_ENUMV(RealSingle), "CoordinateX", xcoor, &Cgns_t::coor_id );
    cg_close( Cgns_t::file_id );
}


void WriteZoneInfo( const std::string & zoneName, ZoneType_t zoneType, cgsize_t * isize )
{
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, zoneName.c_str(), isize, zoneType, & Cgns_t::zone_id );
}


void TestCgnsLink()
{
    if ( cg_open( "testlink.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 1;
    Cgns_t::phys_dim = 3;

    Cgns_t::WriteBase( "Base");

    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    Cgns_t::SetDimArray( dim_array, 10 );
    Cgns_t::SetISize( isize, dim_array );

    ::WriteZoneInfo( "Zone 1", CGNS_ENUMV( Structured ), isize );
    Cgns_t::SetDimArray( dim_array, 20 );
    Cgns_t::SetISize( isize, dim_array );
    ::WriteZoneInfo( "Zone 2", CGNS_ENUMV( Structured ), isize );

    cg_close( Cgns_t::file_id );
}

void cgns_dump_grid()
{
    if ( cg_open( "oneflow-2d.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 2;
    Cgns_t::phys_dim = 3;
    cg_base_write( Cgns_t::file_id, "Base", Cgns_t::cell_dim, Cgns_t::phys_dim, &Cgns_t::base_id );
    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    int ni = 11;
    int nj = 2;
    Cgns_t::SetDimArray( dim_array, ni, nj );
    Cgns_t::SetISize( isize, dim_array );
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, "Zone 1", isize, CGNS_ENUMV(Structured), &Cgns_t::zone_id );
    cgsize_t ipnts[ 6 ];
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1, 1, 1, nj );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "left", CGNS_ENUMV( BCInflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, ni, 1, ni, nj );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "right", CGNS_ENUMV( BCOutflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1, 1, ni, 1 );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "bottom", CGNS_ENUMV( BCSymmetryPlane ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1, nj, ni, nj );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "top", CGNS_ENUMV( BCSymmetryPlane ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }


    cg_close( Cgns_t::file_id );

}

void cgns_dump_grid_1d()
{
    if ( cg_open( "oneflow-1d.cgns", CG_MODE_WRITE, &Cgns_t::file_id ) )
    {
        cg_error_exit();
    }

    Cgns_t::cell_dim = 1;
    Cgns_t::phys_dim = 3;
    cg_base_write( Cgns_t::file_id, "Base", Cgns_t::cell_dim, Cgns_t::phys_dim, &Cgns_t::base_id );
    std::vector<int> dim_array;
    cgsize_t isize[ 9 ];
    int ni = 11;
    Cgns_t::SetDimArray( dim_array, ni );
    Cgns_t::SetISize( isize, dim_array );
    cg_zone_write( Cgns_t::file_id, Cgns_t::base_id, "Zone 1", isize, CGNS_ENUMV(Structured), &Cgns_t::zone_id );
    cgsize_t ipnts[ 6 ];
    {
        Cgns_t::SetStructuredIpnts( ipnts, 1 );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "left", CGNS_ENUMV( BCInflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }
    {
        Cgns_t::SetStructuredIpnts( ipnts, ni );
        cg_boco_write( Cgns_t::file_id, Cgns_t::base_id, Cgns_t::zone_id, "right", CGNS_ENUMV( BCOutflow ),
            CGNS_ENUMV( PointRange ), 2, ipnts, &Cgns_t::bc_id );
    }

    cg_close( Cgns_t::file_id );

}

#endif