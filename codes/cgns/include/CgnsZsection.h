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
#include <vector>
#include <string>
//#include <fstream>
#include "HXCgns.h"
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsSection;
//class ElemFeature;

//class CgnsZone;
//class CgnsZbc;

class CgnsZsection
{
public:
    CgnsZsection( CgnsZone * cgnsZone );
    ~CgnsZsection();
public:
    int nSection;

    HXVector< CgnsSection * > cgnsSections;
    CgnsZone * cgnsZone;
public:
    void AddCgnsSection( CgnsSection * cgnsSection );
    CgnsSection * GetCgnsSection( int iSection );
    void CreateCgnsSection();
    void CreateConnList();
    void ConvertToInnerDataStandard();
    CgnsSection * GetSectionByEid( int eId );
public:
    void ReadNumberOfCgnsSections();
    void ReadCgnsSections();
    void SetElemPosition();
public:
    bool ExistSection( const string & sectionName );
};

#endif

EndNameSpace