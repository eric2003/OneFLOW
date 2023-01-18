/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
    typedef char char33[ 33 ];
public:
    static int file_id, base_id, zone_id, sols_id, field_id;
    static int nbases, nzones, nbccos, nsols, nfields;
    static CGNS_ENUMT( ZoneType_t ) zone_type;
    static int coor_id, bc_id;
    static int cell_dim, phys_dim;
    static cgsize_t isize[ 9 ];
    static cgsize_t irmin[3], irmax[3], cellsize[3];
    static cgsize_t nnodes, ncells;
public:
    static void WriteBase( const std::string & baseName );
    static void SetISize( cgsize_t * isize, std::vector<int> & dim_array );
    static void SetDimArray( std::vector<int> & dim_array, int ni, int nj, int nk );
    static void SetDimArray( std::vector<int> & dim_array, int ni, int nj );
    static void SetDimArray( std::vector<int> & dim_array, int ni );
    static void SetDimensionStr( cgsize_t * isize );
    static void SetStructuredIpnts( cgsize_t * ipnts, int ilo, int jlo, int klo, int ihi, int jhi, int khi );
    static void SetStructuredIpnts( cgsize_t * ipnts, int ilo, int jlo, int ihi, int jhi );
    static void SetStructuredIpnts( cgsize_t * ipnts, int ilo );
};

#endif