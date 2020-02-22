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
#include "HXArray.h"

BeginNameSpace( ONEFLOW )

class BcRecord;

const int F_INNER = 0;
const int F_GHOST = 1;
const int F_TOTAL = 2;

class UnsGrid;
class UGeom
{
public:
    UGeom();
    ~UGeom();
public:
    void Init();
    void CreateBcRegion();
    void SetStEd( int flag );
    void DumpCellFace( int cId );
public:
    UnsGrid * grid;
    int fId;
    int cId;
    int bcfId;
    int nFace;
    int nBFace;
    int nCell;
    int nTCell;
    int ist, ied;
    int lc;
    int rc;
    int ir, bcr, bcdtkey, bctype, nRegion, nRBFace;
    int ireconface;
public:
    IntField * lcf;
    IntField * rcf;
    IntField * blankf;
    LinkField * c2f;

    RealField * xfn;
    RealField * yfn;
    RealField * zfn;
    RealField * vfn;
    RealField * farea;

    RealField * xfc;
    RealField * yfc;
    RealField * zfc;

    RealField * vfx;
    RealField * vfy;
    RealField * vfz;

    RealField * cvol, * cvol1, * cvol2;
    RealField * xcc;
    RealField * ycc;
    RealField * zcc;
    BcRecord * bcRecord;
};

extern UGeom ug;

void AddF2CField( MRField * cellField, MRField * faceField );
void AddF2CFieldDebug( MRField * cellField, MRField * faceField );

EndNameSpace