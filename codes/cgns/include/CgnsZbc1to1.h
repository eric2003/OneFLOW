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

class CgnsZbc1to1
{
public:
    CgnsZbc1to1( CgnsZone * cgnsZone );
    ~CgnsZbc1to1();
public:
    int n1To1;
    HXVector< CgnsBcRegion * > cgnsBcRegion1To1;
    CgnsZone * cgnsZone;
public:
    void AddCgns1To1BcRegion( CgnsBcRegion * cgnsBcRegion );
    CgnsBcRegion * GetCgnsBcRegion1To1( int i1To1 );
    void CreateCgns1To1BcRegion();
    void ConvertToInnerDataStandard();
    void ReadNumberOfCgns1To1();
    void ReadCgns1to1BcRegion();
    void SetPeriodicBc();
};

#endif

EndNameSpace