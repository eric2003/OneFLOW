/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "Grid.h"
#include "Multiarray.h"

BeginNameSpace( ONEFLOW )

class DataBook;
class FaceTopo;
class BcRegionGroup;

class StrGrid : public Grid
{
public:
    IMPLEMENT_GRID_CLONE( StrGrid )
public:
    StrGrid();
    ~StrGrid();
public:
    int  ni, nj, nk;
    FaceTopo * faceTopo;
    BcRegionGroup * bcRegionGroup;
    Field3D * strx, * stry, * strz;
public:
    void Decode( DataBook * databook );
    void Encode( DataBook * databook );
public:
    void ReadGrid( DataBook * databook );
    void WriteGrid( DataBook * databook );
public:
    void ReadBoundaryTopology( DataBook * databook );
    void WriteGridFaceTopology( DataBook * databook );
    void WriteBoundaryTopology( DataBook * databook );
public:
    void SetBasicDimension();
    void SetLayout();
public:
    int ComputeNumberOfNode();
    int ComputeNumberOfCell();
    int ComputeNumberOfFace();
public:
    void GetMinMaxDistance( Real & dismin, Real & dismax );
    void ComputeMinMaxDis3D( Real & dismin, Real & dismax );
    void ComputeMinMaxDis2D( Real & dismin, Real & dismax );
    void ComputeMinMaxDis1D( Real & dismin, Real & dismax );
};

int ComputeNumberOfFace( const int & ni );
int ComputeNumberOfFace( const int & ni, const int & nj );
int ComputeNumberOfFace( const int & ni, const int & nj, const int & nk );

int ComputeNumberOfNode( const int & ni );
int ComputeNumberOfNode( const int & ni, const int & nj );
int ComputeNumberOfNode( const int & ni, const int & nj, const int & nk );

int ComputeNumberOfCell( const int & ni );
int ComputeNumberOfCell( const int & ni, const int & nj );
int ComputeNumberOfCell( const int & ni, const int & nj, const int & nk );

StrGrid * StrGridCast( Grid * gridIn );
void GetIJKRange( int ni, int nj, int nk, int startShift, int endShift, Range & I, Range & J, Range & K );

class IJKRange
{
public:
    IJKRange();
    ~IJKRange();
public:
    static Range I, J, K;
    static int ist, ied;
    static int jst, jed;
    static int kst, ked;
public:
    static void Compute( int ni, int nj, int nk, int ss, int es );
    static void ToScalar();
};

EndNameSpace