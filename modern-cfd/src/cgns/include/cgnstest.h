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

#pragma once
#ifdef PRJ_ENABLE_CGNS
#include <vector>
#include <string>
#include <cgnslib.h>

const int ONE_D = 1;
const int TWO_D = 2;
const int THREE_D = 3;

class CgnsBase;

class Cgns_t
{
public:
    static int file_id, base_id, zone_id;
    static int coor_id, bc_id;
    static int cell_dim, phys_dim;
public:
    static void WriteBase( const std::string & baseName );
    static void SetISize( cgsize_t * isize, std::vector<int> & dim_array );
    static void SetDimArray( std::vector<int> & dim_array, int ni, int nj, int nk );
    static void SetDimArray( std::vector<int> & dim_array, int ni, int nj );
    static void SetDimArray( std::vector<int> & dim_array, int ni );
    static void SetStructuredIpnts( cgsize_t * ipnts, int ilo, int jlo, int klo, int ihi, int jhi, int khi );
    static void SetStructuredIpnts( cgsize_t * ipnts, int ilo, int jlo, int ihi, int jhi );
    static void SetStructuredIpnts( cgsize_t * ipnts, int ilo );
};

class CgnsBase
{
public:
    CgnsBase();
    ~CgnsBase();
public:
    ;
};

void cgns_write_base_test();

void TestCgnsLink();
int cgnstest();
void cgns_dump_grid( float * xcoor, int ni, const std::string & filename );
void cgns_dump_grid();
void cgns_dump_grid_1d();
void fill_isize( cgsize_t * isize, std::vector<int> & dim_array );
#endif