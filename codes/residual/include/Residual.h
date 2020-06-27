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
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

class ResData;

class ResAver
{
public:
    ResAver();
    ~ResAver();
public:
    RealField res;
    int nCell;
public:
    void Init( int nEqu );
    ResAver & operator += ( const ResAver & rhs );
    void Zero();
    void CalcAver( HXVector< ResData > & dataList );
};

class ResMax
{
public:
    ResMax();
    ~ResMax();
public:
    RealField resmax;
    IntField index;
    IntField zid;
    RealField xcc, ycc, zcc, vol;
public:
    void Init( int nEqu );
    void SwapMax( ResMax & rhs );
    void CalcMax( HXVector< ResData > & dataList );
    int CalcMaxId();
};

class ResData
{
public:
    ResData();
    ~ResData();
public:
    ResMax resmax;
    ResAver resave;
    void Init( int nEqu );
};

class Residual
{
public:
    Residual();
    ~Residual();
public:
    virtual void Dump( int sTid ){};
};


EndNameSpace