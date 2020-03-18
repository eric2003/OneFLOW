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

class CgnsZbcBoco
{
public:
    CgnsZbcBoco( CgnsZone * cgnsZone );
    ~CgnsZbcBoco();
public:
    int nBoco;
    HXVector< CgnsBcBoco * > cgnsBcBocos;
    CgnsZone * cgnsZone;
public:
    void AddCgnsBcBoco( CgnsBcBoco * cgnsBcBoco );
    CgnsBcBoco * GetCgnsBc( int iBoco );
    void CreateCgnsZbc();
    void ShiftBcRegion();
    void ConvertToInnerDataStandard();
    void ScanBcFace( FaceSolver * face_solver );
    void PrintZnboco();
    void ReadZnboco();
    void ReadZnboco( int nBoco );
    void ReadCgnsZbcBoco();
    int GetNumberOfActualBcElements();
    void GenerateUnsBcElemConn( CgIntField& bcConn );
};

#endif

EndNameSpace