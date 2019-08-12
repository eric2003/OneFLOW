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
#include "HXDefine.h"
#include "HXCgns.h"

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsCoor
{
public:
    CgnsCoor();
    ~CgnsCoor();
public:
    int ndim;
    IntField nNodeList;
    HXVector< DataType_t > typeList;
    HXVector< void * > coor;
public:
    void * GetCoor( int iCoor ) { return coor[ iCoor ]; };
    void SetAllData( RealField & x, RealField & y, RealField & z );
    void Alloc( int iCoor, int nNode, DataType_t data_type );
public:
    void SetData( int iCoor, DataType_t data_type, Real * var );
    void DeAlloc();
};


#endif

EndNameSpace