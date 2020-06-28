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
#include "Mid.h"
#include "CalcCoor.h"
#include "SimpleDomain.h"
#include <set>
#include <map>
#include <fstream>
using namespace std;
BeginNameSpace( ONEFLOW )

class Block3D;
class BlkMesh;
class SDomain;
class MDomain;
class CalcCoor;
class Face2D;
class MLine;

class MDomain : public DomData
{
public:
    MDomain();
    ~MDomain();
public:
    int pos;
    HXVector< SDomain * > sDomainList;
    CoorMap * coorMap;
public:
    SDomain * FindSDomain( int fid );
    void AddSubDomain( int fid, IntField & lineList, IntField & posList );
    void CalcSubDomainCtrlCoor();
    void CalcCoor();
    int GetNsubDomain();
    void ConstructMultiDomainTopo();
    void ConstructMultiLineToDomainMap();
    void ConstructMultiPointToDomainMap();
    void ConstructMultiPointToPointMap();
    void CreateInpFaceList( HXVector< Face2D * > &facelist );
    void SetBlkBcMesh( Block3D * blk3d );
};

EndNameSpace