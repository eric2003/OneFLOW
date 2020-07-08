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

#include "CgnsFamilyBc.h"
#include "CgnsBase.h"
#include <iostream>
#include <iomanip>
using namespace std;

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
    bcMap = new map< string, int >;
}

void CgnsFamilyBc::Free()
{
    delete bcMap;
}

void CgnsFamilyBc::Register( const string & regionName, int bcType )
{
    map< string, int >::iterator iter = bcMap->find( regionName );
    if ( iter == bcMap->end() )
    {
        ( * CgnsFamilyBc::bcMap )[ regionName ] = bcType;
    }
}

void CgnsFamilyBc::Unregister( const string & regionName )
{
    bcMap->erase( regionName );
}

int CgnsFamilyBc::GetBcType( const string & regionName )
{
    map< string, int >::iterator iter = bcMap->find( regionName );
    if ( iter == bcMap->end() )
    {
        return -1;
    }

    return iter->second;
}

void CgnsFamilyBc::SetFamilyBc( BCType_t & bcType, const string & bcRegionName )
{
    if ( bcType == FamilySpecified )
    {
        int bcTypeFamily = this->GetBcType( bcRegionName );
        bcType = static_cast< BCType_t >( bcTypeFamily );
    }
}

BCType_t CgnsFamilyBc::GetFamilyBcType( const string & bcFamilyName )
{
    int bcTypeOfFamily = this->GetBcType( bcFamilyName );
    BCType_t bcType = static_cast< BCType_t >( bcTypeOfFamily );
    return bcType;
}

void CgnsFamilyBc::ReadFamilySpecifiedBc()
{
    int fileId = cgnsBase->fileId;
    int baseId = cgnsBase->baseId;

    int nFamilies = -1;
    cg_nfamilies( fileId, baseId, & nFamilies );
    cout << "\n";
    cout << "   CGNS nFamilies = " << nFamilies << "\n";
    CgnsTraits::char33 familyName;
    int nBoco = -1;
    int nGeo = -1;
    for ( int iFam = 1; iFam <= nFamilies; ++ iFam )
    {
        cg_family_read( fileId, baseId, iFam, familyName, & nBoco, & nGeo );
        cout << "   iFam = " << iFam;
        cout << " FamilyName = " << setiosflags( ios::left ) << setw( 15 ) << familyName << " nBoco = " << nBoco << " nGeo = " << nGeo << "\n";
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
            cout << "   FamilyBcName = " << setiosflags( ios::left ) << setw( Width ) << familyBcName;
            cout << " CGNS BcType = " << setiosflags( ios::left ) << setw( 5 ) << familyBcType;
            cout << " CGNS BcName = " << setiosflags(ios::left) << setw( stringWidth ) << GetCgnsBcName( familyBcType ) << "\n";
        }
    }

    cout << "\n";
    cout << "\n";
}

#endif

EndNameSpace