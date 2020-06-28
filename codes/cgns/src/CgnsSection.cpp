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

#include "CgnsSection.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "UnitElement.h"
#include "ElementHome.h"
#include "ElemFeature.h"
#include "LogFile.h"

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsSection::CgnsSection( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->connSize = 0;
    this->pos_shift = 0;
}

CgnsSection::~CgnsSection()
{
}

void CgnsSection::ConvertToInnerDataStandard()
{
    this->startId -= 1;
    this->endId   -= 1;

    for ( int iElem = 0; iElem < this->nElement; ++ iElem )
    {
        int e_type = this->eTypeList[ iElem ];
        int npe;
        cg_npe( static_cast< ElementType_t >( e_type ), & npe );
        int pos = ePosList[ iElem ] + this->pos_shift;
        for ( int iNode = 0; iNode < npe; ++ iNode )
        {
            int id = pos + iNode;
            this->connList[ id ] -= 1;
        }
    }
    int kkk = 1;
}

CgInt * CgnsSection::GetAddress( CgInt eId )
{
    int pos = this->ePosList[ eId ] + this->pos_shift;
    return & this->connList[ pos ];
}

void CgnsSection::GetElementNodeId( CgInt eId, CgIntField & eNodeId )
{
    int eNodeNumber  = ONEFLOW::GetElementNodeNumbers( this->eTypeList[ eId ] );
    CgInt * eAddress = this->GetAddress( eId );

    eNodeId.resize( 0 );
    for ( int iNode = 0; iNode < eNodeNumber; ++ iNode )
    {
        eNodeId.push_back( eAddress[ iNode ] );
    }
}

void CgnsSection::SetElementTypeAndNode( ElemFeature * elem_feature )
{
    for ( int iElem = 0; iElem < this->nElement; ++ iElem )
    {
        int e_type = this->eTypeList[ iElem ];

        if ( ! ONEFLOW::IsBasicVolumeElementType( e_type ) ) continue;

        elem_feature->eType->push_back( e_type );

        CgIntField eNodeId;
        this->GetElementNodeId( iElem, eNodeId );

        int eNodeNumber  = ONEFLOW::GetElementNodeNumbers( this->eTypeList[ iElem ] );

        for ( int iNode = 0; iNode < eNodeNumber; ++ iNode )
        {
            eNodeId[ iNode ] = this->cgnsZone->l2g[ eNodeId[ iNode ] ];
        }
        elem_feature->eNodeId.push_back( eNodeId );
    }
}

void CgnsSection::ReadCgnsSection()
{
    this->ReadCgnsSectionInfo();

    this->CreateConnList();

    this->ReadCgnsSectionConnectionList();

    this->SetElemPosition();
}


void CgnsSection::SetSectionInfo( const string & sectionName, int elemType, int startId, int endId )
{
    this->sectionName = sectionName;
    this->eType = elemType;
    this->startId = startId;
    this->endId = endId;
}

void CgnsSection::ReadCgnsSectionInfo()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    ElementType_t elementType;
    CgnsTraits::char33 cgnsSectionName;

    cg_section_read( fileId, baseId, zId, this->id, cgnsSectionName, & elementType, & this->startId, & this->endId, & nbndry, & iparentflag );

    this->sectionName = cgnsSectionName;

    cout << "   Section Name = " << cgnsSectionName << "\n";
    cout << "   Section Type = " << ElementTypeName[ elementType ] << "\n";
    cout << "   startId, endId = " << this->startId << " " << this->endId << "\n";
    this->eType = elementType;

    if ( IsMixedSection() )
    {
        elementDataSize = -1;
        cg_ElementDataSize( fileId, baseId, zId, this->id, & elementDataSize );
        this->pos_shift = 1;
    }
}

void CgnsSection::CreateConnList()
{
    this->CalcNumberOfSectionElements();

    this->CalcCapacityOfCgnsConnectionList();

    this->AllocateCgnsConnectionList();
}

void CgnsSection::CalcNumberOfSectionElements()
{
    this->nElement = this->endId - this->startId + 1;
}

void CgnsSection::CalcCapacityOfCgnsConnectionList()
{
    if ( eType == MIXED ||
         eType == NGON_n ||
         eType == NFACE_n )
    {
        this->connSize = this->elementDataSize;

    }
    else
    {
        UnitElement * unitElement = ElementHome::GetUnitElement( this->eType );
        int nodeNumber = unitElement->GetElementNodeNumbers( this->eType );

        this->connSize = this->nElement * nodeNumber;
    }
}

void CgnsSection::AllocateCgnsConnectionList()
{
    this->connList.resize( this->connSize );
    if ( this->iparentflag )
    {
        this->iparentdata.resize( this->nElement * 4 );
    }
    this->ePosList.resize( this->nElement + 1 );
    this->eTypeList.resize( this->nElement, this->eType );
}

void CgnsSection::ReadCgnsSectionConnectionList()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // Read the connectivity. Again, the node numbering of the 
    // connectivities start at 1. If internally a starting index 
    // of 0 is used ( typical for C-codes ) 1 must be substracted 
    // from the connectivities read. 

    CgInt *addr = NULL;
    if ( this->iparentflag )
    {
        addr = & iparentdata[ 0 ];
    }
    cg_elements_read( fileId, baseId, zId, this->id, & this->connList[ 0 ], addr );
}

void CgnsSection::SetElemPosition()
{
    if ( IsMixedSection() )
    {
        this->SetElemPositionMixed();
    }
    else
    {
        this->SetElemPositionOri();
    }
}

void CgnsSection::SetElemPositionOri()
{
    int pos = 0;
    ePosList[ 0 ] = pos;
    for ( int iElem = 0; iElem < this->nElement; ++ iElem )
    {
        int npe;
        cg_npe( static_cast< ElementType_t >( this->eType ), & npe );
        pos += npe + this->pos_shift;
        ePosList[ iElem + 1 ] = pos;
    }
}

void CgnsSection::SetElemPositionMixed()
{
    //ePosList.resize( this->nElement + 1 );
    //eTypeList.resize( this->nElement );
    int pos = 0;
    ePosList[ 0 ] = pos;
    for ( int iElem = 0; iElem < this->nElement; ++ iElem )
    {
        int e_type = this->connList[ pos ];
        eTypeList[ iElem ] = e_type;
        int npe;
        cg_npe( static_cast< ElementType_t >( e_type ), & npe );
        pos += npe + this->pos_shift;
        ePosList[ iElem + 1 ] = pos;
    }
}

bool CgnsSection::IsMixedSection()
{
    bool flag = ( eType == MIXED || eType == NGON_n || eType == NFACE_n );
    return flag;
}
#endif
EndNameSpace