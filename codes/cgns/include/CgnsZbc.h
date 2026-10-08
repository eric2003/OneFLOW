/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "HXCgns.h"
#include "GridHandles.h"
#include <memory>

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
    explicit CgnsZbc( CgnsZone & cgnsZone );
    ~CgnsZbc();
public:
    std::unique_ptr< CgnsZbcConn > cgnsZbcConn;
    std::unique_ptr< CgnsZbc1to1 > cgnsZbc1to1;
    std::unique_ptr< CgnsZbcBoco > cgnsZbcBoco;

    CgnsZone & cgnsZone;
public:
    CgnsZbcConn & RequireCgnsZbcConn();
    CgnsZbc1to1 & RequireCgnsZbc1to1();
    CgnsZbcBoco & RequireCgnsZbcBoco();
    void ScanBcFace( FaceSolver & faceSolver );
public:
    void ConvertToInnerDataStandard();
    void ReadCgnsGridBoundary();
    void DumpCgnsGridBoundary();

    void FillBcPoints( int * start, int * end, cgsize_t * bcpnts, int dimension );
    void FillBcPoints3D( int * start, int * end, cgsize_t * bcpnts );
    void FillInterface( BcRegion * bcRegion, cgsize_t * ipnts, cgsize_t * ipntsdonor, int * itranfrm, int dimension );
    void FillRegion( TestRegion * r, cgsize_t * ipnts, int dimension );
    void DumpCgnsGridBoundary( Grid * gridIn, const Grids & grids );
public:
    void CreateCgnsZbc( CgnsZbc * cgnsZbcIn );
public:
    void GenerateUnsBcElemConn( CgIntField& bcConn );
    int GetNumberOfActualBcElements();
    void SetPeriodicBc();
};

#endif

EndNameSpace
