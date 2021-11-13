/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#ifdef ENABLE_CGNS

class CgnsBase;
class CgnsZone;
class CgnsFile;

class CgnsZbase
{
public:
    CgnsZbase ();
    ~CgnsZbase();
public:
    CgnsFile * cgnsFile;
    int nBases;
 
    HXVector< CgnsBase * > baseVector;
public:
    int GetSystemZoneType();
    void ReadCgnsGrid( const string & fileName );
    void ReadCgnsMultiBase();
    void DumpCgnsMultiBase();
    void ReadNumCgnsBase();
    void ConvertToInnerDataStandard();
    void ProcessCgnsBases();
    void OpenCgnsFile( const string & fileName, int cgnsOpenMode );
    void CloseCgnsFile();
public:
    void AddCgnsBase( CgnsBase * cgnsBase );
    CgnsBase * CreateCgnsBase();
    void InitCgnsBase();
public:
    int GetNZones();
    CgnsBase * GetCgnsBase( int iBase );
    CgnsZone * GetCgnsZone( int globalZoneId );
    CgnsZone * GetMultiBaseCgnsZone( int iBase, int iZone );
public:
    void FreeCgnsBases();
    void CreateCgnsZones( int nZones = 0 );
    CgnsZone * CreateCgnsZone();
};

#endif

EndNameSpace
