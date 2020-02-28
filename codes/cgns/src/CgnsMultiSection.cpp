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

#include "CgnsMultiSection.h"
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

CgnsMultiSection::CgnsMultiSection( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
}

CgnsMultiSection::~CgnsMultiSection()
{
}

void CgnsMultiSection::AddCgnsSection( CgnsSection * cgnsSection )
{
    this->cgnsSections.push_back( cgnsSection );
    int secId = cgnsSections.size();
    cgnsSection->id = secId;
}

CgnsSection * CgnsMultiSection::GetCgnsSection( int iSection )
{
    return this->cgnsSections[ iSection ];
}

bool CgnsMultiSection::ExistSection( const string & sectionName )
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        if ( cgnsSection->sectionName == sectionName ) return true;
    }
    return false;
}

void CgnsMultiSection::CreateCgnsSection()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = new CgnsSection( cgnsZone );
        this->AddCgnsSection( cgnsSection );
    }
}

void CgnsMultiSection::CreateConnList()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->CreateConnList();
    }
}

void CgnsMultiSection::ConvertToInnerDataStandard()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->ConvertToInnerDataStandard();
    }
}

CgnsSection * CgnsMultiSection::GetSectionByEid( int eId )
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

void CgnsMultiSection::ReadNumberOfCgnsSections()
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

void CgnsMultiSection::ReadCgnsSections()
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

void CgnsMultiSection::SetElemPosition()
{
    for ( int iSection = 0; iSection < this->nSection; ++ iSection )
    {
        CgnsSection * cgnsSection = this->GetCgnsSection( iSection );
        cgnsSection->SetElemPosition();
    }
}

#endif
EndNameSpace