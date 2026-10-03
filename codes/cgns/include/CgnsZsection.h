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
#include <memory>

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class CgnsSection;

class CgnsZsection
{
public:
    explicit CgnsZsection( CgnsZone & cgnsZone );
    ~CgnsZsection();
private:
    HXVector< std::unique_ptr< CgnsSection > > cgnsSections;
    CgnsZone & cgnsZone;
public:
    void AddCgnsSection( std::unique_ptr< CgnsSection > cgnsSection );
    CgnsSection & GetCgnsSection( int iSection );
    const CgnsSection & GetCgnsSection( int iSection ) const;
    int GetNSections() const;
    bool HasPolygonSection() const;
    void CreateCgnsSections( int nSections );
    void CreateConnList();
    void ConvertToInnerDataStandard();
    CgnsSection * GetSectionByEid( int eId );
public:
    int ReadNumberOfCgnsSections();
    void ReadCgnsSections();
    void DumpCgnsSections();
    void SetElemPosition();
public:
    bool ExistSection( const std::string & sectionName );
};

#endif

EndNameSpace
