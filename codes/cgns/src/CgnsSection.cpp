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

#include "CgnsSection.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsFile.h"
#include "StringUtils.h"
#include "Dimension.h"
#include "UnitElement.h"
#include "ElementHome.h"
#include "ElemFeature.h"
#include "LogFile.h"

#include <iostream>
#include <stdexcept>

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsSection::CgnsSection( CgnsZone & cgnsZone )
    : cgnsZone( cgnsZone )
{
    this->connSize = 0;
    this->pos_shift = 0;
    this->nbndry = 0;
    this->iparentflag = 0;
}

CgnsSection::~CgnsSection()
{
}

void CgnsSection::ConvertToInnerDataStandard()
{
    //std::cout << "++++++++++++++++++++++++++++++++++++++++\n";
    //std::cout << "++++++++++++++++++++++++++++++++++++++++\n";
    //std::cout << "ConvertToInnerDataStandard\n";
    //std::cout << "++++++++++++++++++++++++++++++++++++++++\n";
    //std::cout << "++++++++++++++++++++++++++++++++++++++++\n";
    this->startId -= 1;
    this->endId   -= 1;

    if ( this->eType == MIXED )
    {
        for ( int iElem = 0; iElem < this->nElement; ++ iElem )
        {
            int e_type = this->eTypeList[ iElem ];
            int npe;
            cg_npe( static_cast< ElementType_t >( e_type ), & npe );
            int pos = ePosList[ iElem ] + ( this->eType == MIXED ? 1 : this->pos_shift );
            for ( int iNode = 0; iNode < npe; ++ iNode )
            {
                int id = pos + iNode;
                this->connList[ id ] -= 1;
            }
        }
    }
    else if ( this->eType == NFACE_n )
    {
        for ( int i = 0; i < this->connList.size(); ++ i )
        {
            int id = this->connList[ i ];
            int newId = std::abs( id ) - 1;
            if ( id < 0 )
            {
                newId = - std::abs( newId );
            }
            this->connList[ i ] = newId;
        }
    }
    else
    {
        //other cases, include NGON_n
        for ( int i = 0; i < this->connList.size(); ++ i )
        {
            this->connList[ i ] -= 1;
        }
    }


}

CgInt * CgnsSection::GetAddress( CgInt eId )
{
    if ( eId < 0 || eId >= this->nElement ||
         static_cast< size_t >( eId + 1 ) >= this->ePosList.size() ||
         static_cast< size_t >( eId ) >= this->eTypeList.size() )
    {
        throw std::runtime_error( "CgnsSection::GetAddress: element index is out of range" );
    }

    const int eNodeNumber = ONEFLOW::GetElementNodeNumbers( this->eTypeList[ eId ] );
    const CgInt pos = this->ePosList[ eId ] + ( this->eType == MIXED ? 1 : this->pos_shift );
    if ( eNodeNumber <= 0 || pos < 0 ||
         static_cast< size_t >( pos ) > this->connList.size() ||
         static_cast< size_t >( eNodeNumber ) > this->connList.size() - static_cast< size_t >( pos ) )
    {
        throw std::runtime_error( "CgnsSection::GetAddress: element connectivity span is invalid" );
    }

    return & this->connList[ static_cast< size_t >( pos ) ];
}

void CgnsSection::GetElementNodeId( CgInt eId, CgIntField & eNodeId )
{
    if ( eId < 0 || eId >= this->nElement )
    {
        throw std::runtime_error( "CgnsSection::GetElementNodeId: element index is out of range" );
    }

    eNodeId.resize( 0 );
    if ( this->eType != NGON_n )
    {
        if ( static_cast< size_t >( eId ) >= this->eTypeList.size() )
        {
            throw std::runtime_error( "CgnsSection::GetElementNodeId: element type data is incomplete" );
        }

        const int eNodeNumber = ONEFLOW::GetElementNodeNumbers( this->eTypeList[ eId ] );
        CgInt * eAddress = this->GetAddress( eId );

        for ( int iNode = 0; iNode < eNodeNumber; ++ iNode )
        {
            eNodeId.push_back( eAddress[ iNode ] );
        }
    }
    else
    {
        // Polygon offsets delimit a variable-length node list for each face.
        if ( static_cast< size_t >( eId + 1 ) >= this->ePosList.size() )
        {
            throw std::runtime_error( "CgnsSection::GetElementNodeId: NGON offsets are incomplete" );
        }

        const CgInt start = this->ePosList[ eId ];
        const CgInt end = this->ePosList[ eId + 1 ];
        if ( start < 0 || end < start ||
             static_cast< size_t >( end ) > this->connList.size() )
        {
            throw std::runtime_error( "CgnsSection::GetElementNodeId: NGON connectivity span is invalid" );
        }

        for ( CgInt i = start; i < end; ++ i )
        {
            eNodeId.push_back( this->connList[ static_cast< size_t >( i ) ] );
        }
    }
}

void CgnsSection::SetElementTypeAndNode( ElemFeature * elem_feature )
{
    for ( int iElem = 0; iElem < this->nElement; ++ iElem )
    {
        int e_type = this->eTypeList[ iElem ];

        if ( ! ONEFLOW::IsBasicVolumeElementType( e_type ) ) continue;

        elem_feature->eTypes.push_back( e_type );

        CgIntField eNodeId;
        this->GetElementNodeId( iElem, eNodeId );

        const int eNodeNumber = ONEFLOW::GetElementNodeNumbers( this->eTypeList[ iElem ] );
        if ( eNodeId.size() != static_cast< size_t >( eNodeNumber ) )
        {
            throw std::runtime_error( "CgnsSection::SetElementTypeAndNode: element connectivity has an unexpected node count" );
        }

        const auto & localToGlobal = this->cgnsZone.l2g;
        for ( int iNode = 0; iNode < eNodeNumber; ++ iNode )
        {
            const CgInt localNodeId = eNodeId[ iNode ];
            if ( localNodeId < 0 ||
                 static_cast< size_t >( localNodeId ) >= localToGlobal.size() )
            {
                throw std::runtime_error( "CgnsSection::SetElementTypeAndNode: local node index is out of range" );
            }
            eNodeId[ iNode ] = localToGlobal[ static_cast< size_t >( localNodeId ) ];
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

void CgnsSection::DumpCgnsSection()
{
    this->DumpCgnsSectionInfo();

    this->DumpCgnsSectionConnectionList();

}

void CgnsSection::SetSectionInfo( const std::string & sectionName, int elemType, int startId, int endId )
{
    this->sectionName = sectionName;
    this->eType = elemType;
    this->startId = startId;
    this->endId = endId;
}

void CgnsSection::ReadCgnsSectionInfo()
{
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    ElementType_t elementType;
    CgnsTraits::char33 cgnsSectionName;

    cg_section_read( fileId, baseId, zId, this->id, cgnsSectionName, & elementType, & this->startId, & this->endId, & nbndry, & iparentflag );

    this->sectionName = cgnsSectionName;

    std::cout << "   Section Name = " << cgnsSectionName << "\n";
    std::cout << "   Section Type = " << ElementTypeName[ elementType ] << "\n";
    std::cout << "   startId, endId = " << this->startId << " " << this->endId << "\n";
    this->eType = elementType;

    this->elementDataSize = -1;
    cg_ElementDataSize( fileId, baseId, zId, this->id, & this->elementDataSize );

    std::cout << "   elementDataSize = " << elementDataSize << "\n";

    //if ( this->IsMixedSection() )
    //{
    //    this->pos_shift = 1;
    //    //this->conn_offsets.resize(this->nElements+1);
    //    //cg_poly_elements_read ( fileId, baseId, zoneId, sectionId, this->conn.data(), this->conn_offsets.data(), 0 );
    //}
}

void CgnsSection::DumpCgnsSectionInfo()
{
    std::cout << "   Section Name = " << sectionName << "\n";
    std::cout << "   Section Type = " << ElementTypeName[ eType ] << "\n";
    std::cout << "   startId, endId = " << this->startId << " " << this->endId << "\n";
    std::cout << "   nbndry, iparentflag = " << this->nbndry << " " << this->iparentflag << "\n";
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
        UnitElement & unitElement = ElementHome::GetUnitElement( this->eType );
        int nodeNumber = unitElement.GetElementNodeNumbers( this->eType );

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
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    // Read the connectivity. Again, the node numbering of the 
    // connectivities start at 1. If internally a starting index 
    // of 0 is used ( typical for C-codes ) 1 must be substracted 
    // from the connectivities read. 

    CgInt *addr = NULL;
    if ( this->iparentflag )
    {
        addr = & iparentdata[ 0 ];
    }

    cg_elements_read( fileId, baseId, zId, this->id, this->connList.data(), addr );

    if ( this->eType == NGON_n || this->eType == NFACE_n )
    {
        this->pos_shift = 1;
        cg_poly_elements_read ( fileId, baseId, zId, this->id, this->connList.data(), this->ePosList.data(), 0 );
    }
}

void CgnsSection::DumpCgnsSectionConnectionList()
{
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    // write element connectivity
    ElementType_t elementType = static_cast< ElementType_t >( this->eType );
    cg_section_write( fileId, baseId, zId, this->sectionName.c_str(), elementType, this->startId, this->endId, this->nbndry, & this->connList[ 0 ], & this->id );
}

void CgnsSection::SetElemPosition()
{
    if (  this->eType == MIXED )
    {
        this->SetElemPositionMixed();
    }
    else if ( this->eType != NGON_n && this->eType != NFACE_n )
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
    size_t pos = 0;
    ePosList[ 0 ] = 0;
    for ( int iElem = 0; iElem < this->nElement; ++ iElem )
    {
        if ( pos >= this->connList.size() )
        {
            throw std::runtime_error( "CgnsSection::SetElemPositionMixed: missing element type tag" );
        }

        const int e_type = this->connList[ pos ];
        int npe = -1;
        cg_npe( static_cast< ElementType_t >( e_type ), & npe );

        if ( npe <= 0 || static_cast< size_t >( npe ) + 1 > this->connList.size() - pos )
        {
            throw std::runtime_error( "CgnsSection::SetElemPositionMixed: invalid element connectivity span" );
        }

        this->eTypeList[ iElem ] = e_type;
        pos += static_cast< size_t >( npe ) + 1;
        ePosList[ iElem + 1 ] = static_cast< CgInt >( pos );
    }

    if ( pos != this->connList.size() )
    {
        throw std::runtime_error( "CgnsSection::SetElemPositionMixed: unused connectivity entries" );
    }
}

bool CgnsSection::IsMixedSection()
{
    bool flag = ( eType == MIXED || eType == NGON_n || eType == NFACE_n );
    return flag;
}
#endif
EndNameSpace
