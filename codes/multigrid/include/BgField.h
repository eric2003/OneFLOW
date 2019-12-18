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
#include "HXVector.h"
BeginNameSpace( ONEFLOW )

class FieldWrap;
class BasicBgField
{
public:
    BasicBgField ();
    ~BasicBgField();
public:
    typedef HXVector< HXVector< HXVector< FieldWrap * > > > Field3DType;
public:
    Field3DType data;
public:
    void Init();
    void Free();
};

class BgField
{
public:
    BgField ();
    ~BgField();
protected:
    static HXVector< BasicBgField * > data;
public:
    static bool flag;
    static void Init();
    static void Free();
    static FieldWrap * GetFieldWrap( int zid, int sid, int fid, int gl );
};


EndNameSpace