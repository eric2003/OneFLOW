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

class BlkElem
{
public:
    BlkElem();
    ~BlkElem();
public:
    int eType;
    IntField nodeId;
    IntField faceTypeList;
    LinkField faceList;
public:
    void Init( int eType );
    int GetFaceType( int iFace ) const { return faceTypeList[ iFace ]; };
    void PushElementFace( int faceType, int p1, int p2 );
    void PushElementFace( int faceType, int p1, int p2, int p3, int p4 );
};

class BlkElemHome
{
public:
    BlkElemHome();
    ~BlkElemHome();
public:
    HXVector< BlkElem * > elems;
    bool initFlag;
public:
    void Init();
    void Free();
    BlkElem * GetBlkElem( int eType );
};

extern BlkElemHome bbElemHome;


class BlkFace
{
public:
    BlkFace();
    ~BlkFace();
public:
    IntField lblk, rblk;
    IntField lloc, rloc;
};


EndNameSpace