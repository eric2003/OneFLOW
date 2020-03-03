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
#include "HXCgns.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsBase;
class FaceSolver;

class FaceSolver;
class CgnsBcRegion;
class Grid;
class BcRegion;
class TestRegion;

class CgnsZbcConn
{
public:
    CgnsZbcConn( CgnsZone * cgnsZone );
    ~CgnsZbcConn();
public:
    int nConn;
    HXVector< CgnsBcRegion * > cgnsBcRegionConn;
    CgnsZone * cgnsZone;
public:
    void AddCgnsConnBcRegion( CgnsBcRegion * cgnsBcRegion );
    void CreateCgnsConnBcRegion();
};

class CgnsBcRegionProxy
{
public:
    CgnsBcRegionProxy( CgnsZone * cgnsZone );
    ~CgnsBcRegionProxy();
public:
    int nBoco, n1To1, nConn;
    HXVector< CgnsBcRegion * > cgnsBcRegionBoco;
    HXVector< CgnsBcRegion * > cgnsBcRegion1To1;
    HXVector< CgnsBcRegion * > cgnsBcRegionConn;
    CgnsZbcConn * cgnsZbcConn;
    CgnsZone * cgnsZone;
public:
    void ScanBcFace( FaceSolver * face_solver );
public:
    void CreateCgnsBocoBcRegion();
    void CreateCgns1To1BcRegion();
    void CreateCgnsConnBcRegion();

    void ConvertToInnerDataStandard();
    void ShiftBcRegion();
    CgnsBcRegion * GetCgnsBcRegionBoco( int iBoco );
    CgnsBcRegion * GetCgnsBcRegion1To1( int i1To1 );
    CgnsBcRegion * GetCgnsBcRegionConn( int iConn );

    int GetNBocoDynamic();
    void AddCgnsBocoBcRegion( CgnsBcRegion * cgnsBcRegion );
    void AddCgns1To1BcRegion( CgnsBcRegion * cgnsBcRegion );
    void AddCgnsConnBcRegion( CgnsBcRegion * cgnsBcRegion );

    void ReadCgnsGridBoundary();

    void ReadCgnsBocoBcRegion();
    void ReadCgnsConnBcRegion();
    void ReadCgns1to1BcRegion();
    void FillBcPoints( int * start, int * end, cgsize_t * bcpnts, int dimension );
    void FillBcPoints3D( int * start, int * end, cgsize_t * bcpnts );
    void FillInterface( BcRegion * bcRegion, cgsize_t * ipnts, cgsize_t * ipntsdonor, int * itranfrm, int dimension );
    void FillRegion( TestRegion * r, cgsize_t * ipnts, int dimension );
    void DumpCgnsGridBoundary( Grid * gridIn );
public:
    void ReadNumberOfCgnsBoco();
    void ReadNumberOfCgns1To1();
    void ReadNumberOfCgnsConn();
    void CreateCgnsBcRegion( CgnsBcRegionProxy * bcRegionProxyIn );
public:
    void ReconstructStrRegion();
    void GenerateUnsBcElemConn( CgIntField& bcConn );
    int GetNumberOfActualBcElements();
    void SetPeriodicBc();
};

#endif

EndNameSpace