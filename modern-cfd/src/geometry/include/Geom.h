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
#include <vector>

class BoundarySolver;

class Geom_t
{
public:
    Geom_t();
    ~Geom_t();
public:
    static void Init();
public:
    static int ni_ghost;
    static int ni_global;
    static int ni_global_total;
public:
    static std::vector<int> zone_nis;
    static std::vector<int> proc_ids;
    static std::vector<int> zone_ids;
};

class Geom
{
public:
    Geom();
    ~Geom();
public:
    void Init();
    void GenerateGrid();
    void GenerateGrid( int ni, float xmin, float xmax, float * xcoor );
    void ComputeGeom();
public:
    int zoneId;
    int nZones;
    int ni;
    int ni_total;
    float * xcoor_global;
    float * xcoor;
public:
    float xlen;
    float dx;
    float * ds;
    float xmin, xmax;
public:
    BoundarySolver * bcSolver;
};
