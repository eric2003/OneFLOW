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

#include "CgnsZsection.h"
#include "CgnsSection.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "CgnsFile.h"
#include "StringUtils.h"
#include "Dimension.h"
#include "UnitElement.h"
#include "ElementHome.h"
#include "ElemFeature.h"
#include <iostream>
#include <utility>


BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZsection::CgnsZsection( CgnsZone & cgnsZone )
    : cgnsZone( cgnsZone )
{
}

CgnsZsection::~CgnsZsection()
{
}

void CgnsZsection::AddCgnsSection( std::unique_ptr< CgnsSection > cgnsSection )
{
    CgnsSection * section = cgnsSection.get();
    this->cgnsSections.push_back( std::move( cgnsSection ) );
    int secId = cgnsSections.size();
    section->id = secId;
}

CgnsSection * CgnsZsection::GetCgnsSection( int iSection )
{
    return this->cgnsSections[ iSection ].get();
}

const CgnsSection * CgnsZsection::GetCgnsSection( int iSection ) const
{
    return this->cgnsSections[ iSection ].get();
}

int CgnsZsection::GetNSections() const
{
    return static_cast< int >( this->cgnsSections.size() );
}

bool CgnsZsection::ExistSection( const std::string & sectionName )
{
    if ( this->cgnsSections.size() == 0 ) return false;
    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        if ( cgnsSection->sectionName == sectionName ) return true;
    }
    return false;
}

bool CgnsZsection::HasPolygonSection() const
{
    for ( int iSection = 0; iSection < this->cgnsSections.size(); ++ iSection )
    {
        const CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        if ( cgnsSection->eType == NGON_n ||
             cgnsSection->eType == NFACE_n )
        {
            return true;
        }
    }
    return false;
}

void CgnsZsection::CreateCgnsSections( int nSections )
{
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        this->AddCgnsSection( std::make_unique< CgnsSection >( & cgnsZone ) );
    }
}

void CgnsZsection::CreateConnList()
{
    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->CreateConnList();
    }
}

void CgnsZsection::ConvertToInnerDataStandard()
{
    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->ConvertToInnerDataStandard();
    }
}

CgnsSection * CgnsZsection::GetSectionByEid( int eId )
{
    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        if ( cgnsSection->startId <= eId && 
             eId <= cgnsSection->endId )
        {
            return cgnsSection;
        }
    }
    return 0;
}

int CgnsZsection::ReadNumberOfCgnsSections()
{
    int fileId = cgnsZone.cgnsBase->cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase->baseId;
    int zId = cgnsZone.zId;
    int nSections = 0;

    // Determine the number of sections for this zone. Note that
    // surface elements can be stored in a cellVolume zone, but they
    // are NOT taken into account in the number obtained from
    // cg_zone_read.

    cg_nsections( fileId, baseId, zId, & nSections );

    std::cout << "   numberOfCgnsSections = " << nSections << "\n";
    return nSections;
}

void CgnsZsection::ReadCgnsSections()
{
    std::cout << "   Reading Cgns Section Data......\n";
    std::cout << "\n";

    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        std::cout << "-->iSection     = " << iSection << " numberOfCgnsSections = " << nSections << "\n";
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->ReadCgnsSection();
    }
}

void CgnsZsection::DumpCgnsSections()
{
    std::cout << "   Dumping Cgns Section Data......\n";
    std::cout << "\n";

    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        std::cout << "-->iSection     = " << iSection << " numberOfCgnsSections = " << nSections << "\n";
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->DumpCgnsSection();
    }
}

void CgnsZsection::SetElemPosition()
{
    const int nSections = this->GetNSections();
    for ( int iSection = 0; iSection < nSections; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->SetElemPosition();
    }
}

#endif
EndNameSpace
