/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "Configure.h"
#include "HXType.h"
#include "HXDefine.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

class ScalarIFaceIJ
{
public:
    ScalarIFaceIJ() ;
    ~ScalarIFaceIJ();
public:
    int zonei, zonej;
    vector< int > ghostCells;
    vector< int > cells;
};

class GridTopo;
class GridTopos;

class ScalarIFace
{
public:
    ScalarIFace() ;
    ~ScalarIFace();
public:
    vector< ScalarIFaceIJ > data;
    vector<int> ifaces;
    vector<int> zones;
    vector<int> cells;
    int zoneid;
public:
    void GetInterface();
    void CalcInterface( GridTopo * gridTopo );

};

class ScalarIFaces
{
public:
    ScalarIFaces() ;
    ~ScalarIFaces();
public:
    vector< ScalarIFace > data;
public:
    void GetInterface();
    void CalcInterface( GridTopos * gridTopos );

};

EndNameSpace