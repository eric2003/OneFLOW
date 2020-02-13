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

#include "CgnsData.h"
#include "CgnsZone.h"
#include "Dimension.h"
#include "CgnsBase.h"
#include "CgnsBcRegion.h"
#include "CgnsBcRegionProxy.h"
#include "StrUtil.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsData::CgnsData()
{
}

CgnsData::~CgnsData()
{
}

void CgnsData::Create( int nSection )
{
    this->nSection = nSection;
    this->startId.resize( nSection );
    this->endId.resize( nSection );
    this->elemType.resize( nSection );
    this->sectionNameList.resize( nSection );
}

void CgnsData::SetDefaultSectionName()
{
    for ( int iSection = 0; iSection < nSection; ++ iSection )
    {
        string sectionName = AddString( "Section", iSection + 1 );
        this->sectionNameList[ iSection ] = sectionName;
    }


}

void CgnsData::FillCgnsData( CgnsZone * cgnsZone )
{
    int nSection = 2;

    this->Create( nSection );
    this->SetDefaultSectionName();

    CgIntField & startId = this->startId;
    CgIntField & endId = this->endId;
    IntField & elemType = this->elemType;

    int nBFace        = 0;
    int nActualBcFace = 0;

    int nBcRegion = cgnsZone->bcRegionProxy->nBcRegion;

    for ( int iBcRegion = 0; iBcRegion < nBcRegion; ++ iBcRegion )
    {
        CgnsBcRegion * cgnsBcRegion = cgnsZone->bcRegionProxy->GetBcRegion( iBcRegion );
        int nBcElement       = cgnsBcRegion->nElements;
        int nActualBcElement = cgnsBcRegion->GetActualNumberOfBoundaryElements();
        nBFace        += nBcElement;
        nActualBcFace += nActualBcElement;

        cout << " iBcRegion  = " << iBcRegion << " numberOfBoundaryElements       = " << nBcElement << "\n";
        cout << " iBcRegion  = " << iBcRegion << " numberOfActualBoundaryElements = " << nActualBcElement << "\n";
    }
    cout << " numberOfBoundaryFaces       = " << nBFace       << "\n";
    cout << " numberOfActualBoundaryFaces = " << nActualBcFace << "\n";
    
    startId[ 0 ] = 1;
    endId[ 0 ] = cgnsZone->nCell;

    startId[ 1 ] = cgnsZone->nCell + 1;
    endId  [ 1 ] = cgnsZone->nCell + nActualBcFace;

    int celldim = cgnsZone->cgnsBase->celldim;

    if ( celldim == ONE_D )
    {
        elemType[ 0 ]  = CGNS_ENUMV( BAR_2 );
        elemType[ 1 ]  = CGNS_ENUMV( NODE );
    }
    else if ( celldim == TWO_D )
    {
        elemType[ 0 ]  = CGNS_ENUMV( QUAD_4 );
        elemType[ 1 ]  = CGNS_ENUMV( BAR_2  );
    }
    else if ( celldim == THREE_D )
    {
        elemType[ 0 ]  = CGNS_ENUMV( HEXA_8 );
        elemType[ 1 ]  = CGNS_ENUMV( QUAD_4 );
    }
}

#endif
EndNameSpace