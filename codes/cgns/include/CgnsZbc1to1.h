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
#include "HXDefine.h"
#include "HXCgns.h"
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
class CgnsBc1to1;

class CgnsZbc1to1
{
public:
    explicit CgnsZbc1to1( CgnsZone & cgnsZone );
    ~CgnsZbc1to1();
private:
    int n1to1;
    HXVector< std::unique_ptr< CgnsBc1to1 > > cgnsBc1to1s;

public:
    int GetN1to1() const;
    CgnsZone & cgnsZone;
public:
    void AddCgns1To1BcRegion( std::unique_ptr< CgnsBc1to1 > cgnsBc1to1 );
    CgnsBc1to1 & GetCgnsBcRegion1to1( int i1to1 );
    void CreateCgnsZbc();
    void ConvertToInnerDataStandard();
    void ReadZn1to1( int n1to1 );
    void ReadZn1to1();
    void PrintZn1to1();
    void ReadCgnsZbc1to1();
    void DumpCgnsZbc1to1();
};

#endif

EndNameSpace
