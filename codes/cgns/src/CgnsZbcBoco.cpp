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

#include "CgnsZbcBoco.h"
#include "CgnsBcBoco.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsFile.h"
#include "Boundary.h"
#include "StringUtils.h"
#include "Dimension.h"
#include "HXMath.h"
#include "HXStd.h"
#include "StrRegion.h"
#include "StrGrid.h"
#include "GridMediator.h"
#include "FaceSolver.h"
#include "BcRecord.h"
#include <iostream>
#include <utility>



BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZbcBoco::CgnsZbcBoco( CgnsZone & cgnsZone )
    : cgnsZone( cgnsZone )
{
    this->nBoco = 0;
}

CgnsZbcBoco::~CgnsZbcBoco() = default;

int CgnsZbcBoco::GetNBoco() const
{
    return this->nBoco;
}

void CgnsZbcBoco::AddCgnsBcBoco( CgnsBcBoco * cgnsBcBoco )
{
    this->AddCgnsBcBoco( std::unique_ptr< CgnsBcBoco >( cgnsBcBoco ) );
}

void CgnsZbcBoco::AddCgnsBcBoco( std::unique_ptr< CgnsBcBoco > cgnsBcBoco )
{
    CgnsBcBoco * bcBoco = cgnsBcBoco.get();
    this->cgnsBcBocos.push_back( std::move( cgnsBcBoco ) );
    int id = this->cgnsBcBocos.size();
    bcBoco->bcId = id;
}

CgnsBcBoco * CgnsZbcBoco::GetCgnsBc( int iBoco )
{
    return this->cgnsBcBocos[ iBoco ].get();
}

void CgnsZbcBoco::CreateCgnsZbc()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
this->AddCgnsBcBoco( std::make_unique< CgnsBcBoco >( &this->cgnsZone ) );
    }
}

void CgnsZbcBoco::ShiftBcRegion()
{
    int baseFlag = 1;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        if ( ! cgnsBcBoco->CalcBase() )
        {
            baseFlag = 0;
            break;
        }
    }

    if ( baseFlag == 0 )
    {
        for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
        {
            CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
            cgnsBcBoco->ShiftBcRegion();
        }
    }
}

void CgnsZbcBoco::ConvertToInnerDataStandard()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        cgnsBcBoco->ConvertToInnerDataStandard();
    }
}

void CgnsZbcBoco::ScanBcFace( FaceSolver & faceSolver )
{
    std::cout << " Now ScanBcFace......\n\n";
    std::cout << " nBoco = " << this->nBoco << std::endl;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        std::cout << " iBoco = " << iBoco << " ";
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        std::cout << " BCTypeName = " << ONEFLOW::GetCgnsBcName( cgnsBcBoco->bcType ) << std::endl;
        std::cout << " BCRegion Name = " << cgnsBcBoco->name << std::endl;

        RegionNameMap::AddRegion( cgnsBcBoco->name );
        int bcNameId = RegionNameMap::FindRegionId( cgnsBcBoco->name );
        cgnsBcBoco->nameId = bcNameId;
        cgnsBcBoco->ScanBcFace( faceSolver );
    }
}

void CgnsZbcBoco::PrintZnboco()
{
    std::cout << "   nBoco        = " << this->nBoco << std::endl;
}

void CgnsZbcBoco::ReadZnboco()
{
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    // Determine the number of boundary conditions for this zone.
    cg_nbocos( fileId, baseId, zId, & this->nBoco );
    this->PrintZnboco();
}

void CgnsZbcBoco::ReadZnboco( int nBoco )
{
    this->nBoco = nBoco;
    this->PrintZnboco();
}

void CgnsZbcBoco::ReadCgnsZbcBoco()
{
    this->ReadZnboco();
    this->CreateCgnsZbc();

    for ( int iBoco = 0; iBoco < nBoco; ++ iBoco )
    {
        std::cout << "\n";
        std::cout << "-->iBoco  = " << iBoco << " nBoco = " << nBoco << "\n";
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        cgnsBcBoco->ReadCgnsBcBoco();
    }
}

void CgnsZbcBoco::DumpCgnsZbcBoco()
{
    this->PrintZnboco();

    for ( int iBoco = 0; iBoco < nBoco; ++ iBoco )
    {
        std::cout << "\n";
        std::cout << "-->iBoco  = " << iBoco << " nBoco = " << nBoco << "\n";
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        cgnsBcBoco->DumpCgnsBcBoco();
    }
}

CgnsBcBoco * CgnsZbcBoco::WriteCgnsBoco( const std::string & bocoName, BCType_t bocotype,  PointSetType_t ptset_type, cgsize_t npnts, const cgsize_t * pnts )
{
    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zId = cgnsZone.zId;

    std::unique_ptr< CgnsBcBoco > ownedBcBoco = std::make_unique< CgnsBcBoco >( &this->cgnsZone );
    CgnsBcBoco * cgnsBcBoco = ownedBcBoco.get();
    this->AddCgnsBcBoco( std::move( ownedBcBoco ) );

    cgnsBcBoco->WriteCgnsBoco( bocoName, bocotype, ptset_type, npnts, pnts );

    return cgnsBcBoco;
}

int CgnsZbcBoco::GetNumberOfActualBcElements()
{
    int nBFaces = 0;
    int nActualBcFace = 0;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        int nBcElement = cgnsBcBoco->nElements;
        int nActualBcElement = cgnsBcBoco->GetActualNumberOfBoundaryElements();
        nBFaces += nBcElement;
        nActualBcFace += nActualBcElement;

        std::cout << " iBoco  = " << iBoco << " numberOfBoundaryElements       = " << nBcElement << "\n";
        std::cout << " iBoco  = " << iBoco << " numberOfActualBoundaryElements = " << nActualBcElement << "\n";
    }
    std::cout << " numberOfBoundaryFaces       = " << nBFaces << "\n";
    std::cout << " numberOfActualBoundaryFaces = " << nActualBcFace << "\n";
    return nActualBcFace;
}

void CgnsZbcBoco::GenerateUnsBcElemConn( CgIntField& bcConn )
{
    int nBcElem = 0;
    int pos = 0;

    std::cout << " pos = " << pos << "\n";

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * bcRegion = this->GetCgnsBc( iBoco );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );
        SetBcConn( &cgnsZone, ijkMin, ijkMax, bcConn, pos, nBcElem );
        std::cout << " pos = " << pos << "\n";
        std::cout << " nBcElem = " << nBcElem << " boundaryElementSize = " << nBcElem * 4 << "\n";
    }
}


#endif
EndNameSpace
