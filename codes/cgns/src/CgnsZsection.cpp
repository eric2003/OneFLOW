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

#include "CgnsZsection.h"
#include "CgnsSection.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "UnitElement.h"
#include "ElementHome.h"
#include "ElemFeature.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZsection::CgnsZsection( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
}

CgnsZsection::~CgnsZsection()
{
}

void CgnsZsection::AddCgnsSection( CgnsSection * cgnsSection )
{
    this->cgnsSections.push_back( cgnsSection );
    int secId = cgnsSections.size();
    cgnsSection->id = secId;
}

CgnsSection * CgnsZsection::GetCgnsSection( int iSection )
{
    return this->cgnsSections[ iSection ];
}

bool CgnsZsection::ExistSection( const string & sectionName )
{
    if ( this->cgnsSections.size() == 0 ) return false;
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        if ( cgnsSection->sectionName == sectionName ) return true;
    }
    return false;
}

void CgnsZsection::CreateCgnsSection()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = new CgnsSection( cgnsZone );
        this->AddCgnsSection( cgnsSection );
    }
}

void CgnsZsection::CreateConnList()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->CreateConnList();
    }
}

void CgnsZsection::ConvertToInnerDataStandard()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->ConvertToInnerDataStandard();
    }
}

CgnsSection * CgnsZsection::GetSectionByEid( int eId )
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
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

void CgnsZsection::ReadNumberOfCgnsSections()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // Determine the number of sections for this zone. Note that
    // surface elements can be stored in a cellVolume zone, but they
    // are NOT taken into account in the number obtained from 
    // cg_zone_read.

    cg_nsections( fileId, baseId, zId, & this->nSection );

    cout << "   numberOfCgnsSections = " << this->nSection << "\n";
}

void CgnsZsection::ReadCgnsSections()
{
    cout << "   Reading Cgns Section Data......\n";
    cout << "\n";

    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        cout << "-->iSection     = " << iSection << " numberOfCgnsSections = " << this->nSection << "\n";
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->ReadCgnsSection();
    }
}

void CgnsZsection::SetElemPosition()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->SetElemPosition();
    }
}

#endif
EndNameSpace