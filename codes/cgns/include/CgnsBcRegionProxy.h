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
#include <vector>
#include <string>
#include <fstream>
#include "HXCgns.h"
using namespace std;

BeginNameSpace( ONEFLOW )

class CgnsZone;
class CgnsBase;
class FaceSolver;

class FaceSolver;
class CgnsBcRegion;

class CgnsBcRegionProxy
{
public:
	CgnsBcRegionProxy( CgnsZone * cgnsZone );
	~CgnsBcRegionProxy();
public:
	int nBcRegion, n1To1, nOrdinaryBcRegion;
	int n1To1General;
	int nConn;
	HXVector< CgnsBcRegion * > cgnsBcRegions;
	HXVector< CgnsBcRegion * > bcRegion1To1;
	CgnsZone * cgnsZone;
public:
    void ScanBcFace( FaceSolver * face_solver );
public:
    void CreateCgnsBcRegion();
    void ConvertToInnerDataStandard();
    CgnsBcRegion * GetBcRegion( int ir );

    void ReadCgnsOrdinaryBcRegion();
    void ReadCgnsGridBoundary();
    void ReadCgnsInterfaceBcRegion();
public:
    void ReadNumberCgnsConnBcInfo();
    void ReadNumberOfCgnsOrdinaryBcRegions();
    void ReadNumberOfCgns1To1BcRegions();
    void ReadNumberOfCgnsConn();
    void CreateCgnsBcRegion( CgnsBcRegionProxy * bcRegionProxyIn );
public:
	void ReconstructStrRegion();
    void GenerateUnsBcElemConn( IntField & bcConn );
	void SetPeriodicBc();
};


EndNameSpace