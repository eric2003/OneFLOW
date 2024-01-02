/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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

#include "CgnsFamilyBc.h"
#include "CgnsBase.h"
#include "CgnsFile.h"
#include <iostream>
#include <iomanip>

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsFamilyBc::CgnsFamilyBc( CgnsBase * cgnsBase )
{
    this->cgnsBase = cgnsBase;
    Init();
}

CgnsFamilyBc::~CgnsFamilyBc()
{
    Free();
}

void CgnsFamilyBc::Init()
{
    bcMap = new std::map< std::string, int >;
}

void CgnsFamilyBc::Free()
{
    delete bcMap;
}

void CgnsFamilyBc::Register( const std::string & regionName, int bcType )
{
    std::map< std::string, int >::iterator iter = bcMap->find( regionName );
    if ( iter == bcMap->end() )
    {
        ( * CgnsFamilyBc::bcMap )[ regionName ] = bcType;
    }
}

void CgnsFamilyBc::Unregister( const std::string & regionName )
{
    bcMap->erase( regionName );
}

int CgnsFamilyBc::GetBcType( const std::string & regionName )
{
    std::map< std::string, int >::iterator iter = bcMap->find( regionName );
    if ( iter == bcMap->end() )
    {
        return -1;
    }

    return iter->second;
}

void CgnsFamilyBc::SetFamilyBc( BCType_t & bcType, const std::string & bcRegionName )
{
    if ( bcType == FamilySpecified )
    {
        int bcTypeFamily = this->GetBcType( bcRegionName );
        bcType = static_cast< BCType_t >( bcTypeFamily );
    }
}

BCType_t CgnsFamilyBc::GetFamilyBcType( const std::string & bcFamilyName )
{
    int bcTypeOfFamily = this->GetBcType( bcFamilyName );
    BCType_t bcType = static_cast< BCType_t >( bcTypeOfFamily );
    return bcType;
}

void CgnsFamilyBc::ReadFamilySpecifiedBc()
{
    int fileId = cgnsBase->cgnsFile->fileId;
    int baseId = cgnsBase->baseId;

    int nFamilies = -1;
    cg_nfamilies( fileId, baseId, & nFamilies );
    std::cout << "\n";
    std::cout << "   CGNS nFamilies = " << nFamilies << "\n";
    CgnsTraits::char33 familyName;
    int nBoco = -1;
    int nGeo = -1;
    for ( int iFam = 1; iFam <= nFamilies; ++ iFam )
    {
        cg_family_read( fileId, baseId, iFam, familyName, & nBoco, & nGeo );
        std::cout << "   iFam = " << iFam;
        std::cout << " FamilyName = " << std::setiosflags( std::ios::left ) << std::setw( 15 ) << familyName << " nBoco = " << nBoco << " nGeo = " << nGeo << "\n";
    }

    for ( int iFam = 1; iFam <= nFamilies; ++ iFam )
    {
        cg_family_read( fileId, baseId, iFam, familyName, & nBoco, & nGeo );
        if ( nBoco == 1 )
        {
            CgnsTraits::char33 familyBcName;
            BCType_t familyBcType = BCTypeNull;
            cg_fambc_read( fileId, baseId, iFam, nBoco, familyBcName, & familyBcType );
            this->Register( familyName, familyBcType );
            int Width = 10;
            int stringWidth = 23;
            std::cout << "   FamilyBcName = " << std::setiosflags( std::ios::left ) << std::setw( Width ) << familyBcName;
            std::cout << " CGNS BcType = " << std::setiosflags( std::ios::left ) << std::setw( 5 ) << familyBcType;
            std::cout << " CGNS BcName = " << std::setiosflags(std::ios::left) << std::setw( stringWidth ) << GetCgnsBcName( familyBcType ) << "\n";
        }
    }

    std::cout << "\n";
    std::cout << "\n";
}

#endif

EndNameSpace
