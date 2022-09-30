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

#include "Cgns_t.h"
#include "cgnslib.h"
#include "tools.h"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#ifdef PRJ_ENABLE_CGNS

int Cgns_t::cell_dim = 3;
int Cgns_t::phys_dim = 3;
int Cgns_t::file_id = -1;
int Cgns_t::base_id = -1;
int Cgns_t::zone_id = -1;
int Cgns_t::coor_id = -1;
int Cgns_t::bc_id   = -1;
int Cgns_t::nbases  = -1;
int Cgns_t::nzones  = -1;
int Cgns_t::nbccos  = -1;
int Cgns_t::nsols   = -1;
int Cgns_t::sols_id = -1;
int Cgns_t::nfields = -1;
int Cgns_t::field_id = -1;

CGNS_ENUMT( ZoneType_t ) Cgns_t::zone_type = CGNS_ENUMV( ZoneTypeNull );
cgsize_t Cgns_t::isize[ 9 ];
cgsize_t Cgns_t::irmin[ 3 ];
cgsize_t Cgns_t::irmax[ 3 ];
cgsize_t Cgns_t::cellsize[ 3 ];
cgsize_t Cgns_t::nnodes = -1;
cgsize_t Cgns_t::ncells = -1;

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

void Cgns_t::SetDimensionStr( cgsize_t * isize )
{
    // lower range index
    Cgns_t::irmin[ 0 ] = 1;
    Cgns_t::irmin[ 1 ] = 1;
    Cgns_t::irmin[ 2 ] = 1;

    // upper range index of vertices
    Cgns_t::irmax[ 0 ] = 1;
    Cgns_t::irmax[ 1 ] = 1;
    Cgns_t::irmax[ 2 ] = 1;

    Cgns_t::cellsize[ 0 ] = 1;
    Cgns_t::cellsize[ 1 ] = 1;
    Cgns_t::cellsize[ 2 ] = 1;

    // upper range index of vertices
    // vertex size
    int j = 0;
    for ( int idim = 0; idim < Cgns_t::cell_dim; ++ idim )
    {
        Cgns_t::irmax[ idim ] = isize[ j ++ ];
    }
    // cell size
    for ( int idim = 0; idim < Cgns_t::cell_dim; ++ idim )
    {
        Cgns_t::cellsize[ idim ] = isize[ j ++ ];
    }

    std::cout << "   The Dimension Of Grid is : \n";
    std::cout << "   I Direction " << std::setw( 10 ) << irmin[ 0 ] << std::setw( 10 ) << irmax[ 0 ] << "\n";
    std::cout << "   J Direction " << std::setw( 10 ) << irmin[ 1 ] << std::setw( 10 ) << irmax[ 1 ] << "\n";
    std::cout << "   K Direction " << std::setw( 10 ) << irmin[ 2 ] << std::setw( 10 ) << irmax[ 2 ] << "\n";
    Cgns_t::nnodes = Cgns_t::irmax[ 0 ] * Cgns_t::irmax[ 1 ] * Cgns_t::irmax[ 2 ];
    Cgns_t::ncells = Cgns_t::cellsize[ 0 ] * Cgns_t::cellsize[ 1 ] * Cgns_t::cellsize[ 2 ];
    std::cout << "   Cgns_t::nnodes =  " << Cgns_t::nnodes << "   Cgns_t::ncells =  "  << Cgns_t::ncells << "\n";
}

#endif