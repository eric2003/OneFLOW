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
class CgnsBcBoco;
class Grid;
class BcRegion;
class TestRegion;

class CgnsZbcConn;
class CgnsZbc1to1;
class CgnsZbcBoco;

class CgnsZbc
{
public:
    CgnsZbc( CgnsZone * cgnsZone );
    ~CgnsZbc();
public:
    CgnsZbcConn * cgnsZbcConn;
    CgnsZbc1to1 * cgnsZbc1to1;
    CgnsZbcBoco * cgnsZbcBoco;

    CgnsZone * cgnsZone;
public:
    void ScanBcFace( FaceSolver * face_solver );
public:
    void ConvertToInnerDataStandard();
    void ReadCgnsGridBoundary();

    void FillBcPoints( int * start, int * end, cgsize_t * bcpnts, int dimension );
    void FillBcPoints3D( int * start, int * end, cgsize_t * bcpnts );
    void FillInterface( BcRegion * bcRegion, cgsize_t * ipnts, cgsize_t * ipntsdonor, int * itranfrm, int dimension );
    void FillRegion( TestRegion * r, cgsize_t * ipnts, int dimension );
    void DumpCgnsGridBoundary( Grid * gridIn );
public:
    void CreateCgnsZbc( CgnsZbc * cgnsZbcIn );
public:
    void GenerateUnsBcElemConn( CgIntField& bcConn );
    int GetNumberOfActualBcElements();
    void SetPeriodicBc();
};

#endif

EndNameSpace