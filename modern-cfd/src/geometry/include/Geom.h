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
#include <string>

class BoundarySolver;

class Geom_t
{
public:
    Geom_t();
    ~Geom_t();
public:
    static void Init();
    static void Finalize();
    static void DumpGrid( const std::string & fileName );
    static void GenerateGrid();
    static void ReadGrid( const std::string & gridName );
public:
    static int ni_ghost;
    static int ni_global;
    static int ni_global_total;
    static float * xcoor_global;
    static float dx;
public:
    static std::vector<int> zone_nis;
    static std::vector<int> proc_ids;
    static std::vector<int> zone_ids;
public:
};

class Geom
{
public:
    Geom();
    ~Geom();
public:
    void Init();
    void GenerateGrid();
    void ComputeGeom();
public:
    int zoneId;
    int nZones;
    int ni;
    int ni_total;
    float * xcoor;
public:
    float * ds;
public:
    BoundarySolver * bcSolver;
};
