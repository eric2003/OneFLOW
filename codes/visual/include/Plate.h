/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "HXClone.h"
#include "Task.h"
#include "HXArray.h"
#include <sstream>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

const int X_DIR = 0;
const int Y_DIR = 1;
const int Z_DIR = 2;

class SliceInfo
{
public:
    SliceInfo();
    ~SliceInfo();
public:
    RealField slicepos;
    IntField dir1, dir2;
public:
    void AddSlice( Real pos, int dir1, int dir2 );
};

class PlaneData
{
public:
    PlaneData();
    ~PlaneData();
public:
    MRField * nodedata;
    RealField2D slicedata;
    RealField x, y, z;
    RealField var1, var2, var;
public:
    void Init();
    int  GetNNode() { return static_cast<int> (x.size()); }
public:
    void AddPoint( Real xm, Real ym, Real zm );
    void AddVar( RealField & var );
    void AddVar( int p1, int p2, Real c1, Real c2 );
    void Write( DataBook * dataBook );
    void Read( DataBook * dataBook );
    void SortData( RealField & varList );
    void SortDataByAxis( int axis );
    int FindYIndex();
};

class LamData
{
public:
    LamData();
    ~LamData();
public:
    HXVector< PlaneData * > data;
public:
    void Init();
    void AddVar( RealField & point, int p1, int p2, Real c1, Real c2 );
    int  GetNNode();
    void Write( DataBook * dataBook );
    void Read( DataBook * dataBook );
    void SortDataByAxis( int axis );
    int FindYIndex();
};

class UnsGrid;

class CuttingClass
{
public:
    CuttingClass();
    ~CuttingClass();
public:
    StringField nameList;
    SliceInfo sliceInfo;
    HXVector< LamData * > sliceData;
public:
    void Init();
    void Slice();
    void Swap();
    virtual void Dump(){};
    virtual void Dump( LamData * lamData, fstream & file ){};
    void Write( DataBook * dataBook );
    void Read( DataBook * dataBook );
    RealField & GetCoor( UnsGrid * grid, int cutAxis );
    void CutPlane( Real cutPosition, int cutAxis, LamData * lamData );
};


EndNameSpace