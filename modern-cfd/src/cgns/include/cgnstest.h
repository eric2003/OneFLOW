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

class CgnsBase;

class CgnsZone
{
public:
    CgnsZone() {};
    ~CgnsZone() {};
public:
    cgsize_t isize[ 9 ];
    cgsize_t irmin[3], irmax[3];
};

void cgns_write_base_test();

void TestCgnsLink();
int cgnstest();
void ReadFieldExample();
void ModifyFieldExample();
void WriteFieldExample( int ni );
void cgns_read_grid( std::vector<float> & xcoor, const std::string & filename );
void cgns_read_and_modify_grid( std::vector<float> & xcoor, const std::string & filename );
void cgns_dump_grid( float * xcoor, int ni, const std::string & filename );
void cgns_dump_grid();
void cgns_dump_grid_1d();
#endif