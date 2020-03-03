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

#include "CgnsZbcBoco.h"
#include "CgnsBcBoco.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "Boundary.h"
#include "StrUtil.h"
#include "Dimension.h"
#include "HXMath.h"
#include "HXStd.h"
#include "StrRegion.h"
#include "StrGrid.h"
#include "GridMediator.h"
#include "FaceSolver.h"
#include "BcRecord.h"
#include <iostream>

using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZbcBoco::CgnsZbcBoco( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->nBoco = 0;
}

CgnsZbcBoco::~CgnsZbcBoco()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        delete this->cgnsBcRegionBoco[ iBoco ];
    }
}

void CgnsZbcBoco::AddCgnsBocoBcRegion( CgnsBcBoco * cgnsBcRegion )
{
    this->cgnsBcRegionBoco.push_back( cgnsBcRegion );
}

CgnsBcBoco * CgnsZbcBoco::GetCgnsBcRegionBoco( int iBoco )
{
    return this->cgnsBcRegionBoco[ iBoco ];
}

void CgnsZbcBoco::CreateCgnsBocoBcRegion()
{
    cout << "   nBoco        = " << this->nBoco << endl;
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcRegion = new CgnsBcBoco( this->cgnsZone );
        this->AddCgnsBocoBcRegion( cgnsBcRegion );
    }
}

void CgnsZbcBoco::ShiftBcRegion()
{
    int baseFlag = 1;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
        if ( ! cgnsBcRegion->ComputeBase() )
        {
            baseFlag = 0;
            break;
        }
    }

    if ( baseFlag == 0 )
    {
        for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
        {
            CgnsBcBoco * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
            cgnsBcRegion->ShiftBcRegion();
        }
    }
}

void CgnsZbcBoco::ConvertToInnerDataStandard()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
        cgnsBcRegion->ConvertToInnerDataStandard();
    }
}

void CgnsZbcBoco::ScanBcFace( FaceSolver * face_solver )
{
    cout << " Now ScanBcFace......\n\n";
    cout << " nBoco = " << this->nBoco << endl;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        cout << " iBoco = " << iBoco << " ";
        CgnsBcBoco * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
        cout << " BCTypeName = " << ONEFLOW::GetCgnsBcName( cgnsBcRegion->bcType ) << endl;
        cout << " BCRegion Name = " << cgnsBcRegion->name << endl;

        RegionNameMap::AddRegion( cgnsBcRegion->name );
        int bcNameId = RegionNameMap::FindRegionId( cgnsBcRegion->name );
        cgnsBcRegion->nameId = bcNameId;
        cgnsBcRegion->bcId = iBoco + 1;
        cgnsBcRegion->ScanBcFace( face_solver );
    }
    face_solver->ScanInterfaceBc();
}

void CgnsZbcBoco::ReadNumberOfCgnsBoco()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // Determine the number of boundary conditions for this zone.
    cg_nbocos( fileId, baseId, zId, & this->nBoco );
}

void CgnsZbcBoco::ReadCgnsBocoBcRegion()
{
    this->ReadNumberOfCgnsBoco();
    this->CreateCgnsBocoBcRegion();

    for ( int iBcRegion = 0; iBcRegion < nBoco; ++ iBcRegion )
    {
        cout << "\n-->iBcRegion  = " << iBcRegion;
        cout << " nOrdinaryBcRegion = " << nBoco << "\n";
        CgnsBcBoco * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBcRegion );
        cgnsBcRegion->bcId = iBcRegion + 1;
        cgnsBcRegion->ReadCgnsBocoBcRegion();
    }
}

void CgnsZbcBoco::ReconstructStrRegion()
{
    int ni = static_cast<int> (this->cgnsZone->GetNI());
    int nj = static_cast<int> (this->cgnsZone->GetNJ());
    int nk = static_cast<int> (this->cgnsZone->GetNK());

    if ( nk == 1 ) return;

    MyRegionFactory rfact;
    rfact.ni = ni;
    rfact.nj = nj;
    rfact.nk = nk;
    rfact.CreateRegion();

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * bcRegion = this->GetCgnsBcRegionBoco( iBoco );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );

        rfact.AddRefBcRegion( ijkMin, ijkMax );
    }

    rfact.Run();

    UInt nnr = rfact.bcregions.size();

    nBoco += static_cast<int> (nnr);

    for ( UInt i = 0; i < nnr; ++ i )
    {
        CgnsBcBoco * rr = new CgnsBcBoco( this->cgnsZone );
        MyRegion * r = rfact.bcregions[ i ];

        int id = static_cast<int> ( this->GetNBocoDynamic() + 1 );
        rr->bcId = id;
        rr->ReconstructStrRegion( r->ijkmin, r->ijkmax );

        this->AddCgnsBocoBcRegion( rr );
    }
    int kkk = 1;
}

int CgnsZbcBoco::GetNBocoDynamic()
{
    return this->cgnsBcRegionBoco.size();
}

int CgnsZbcBoco::GetNumberOfActualBcElements()
{
    int nBFace = 0;
    int nActualBcFace = 0;

    for ( int iBcRegion = 0; iBcRegion < this->nBoco; ++ iBcRegion )
    {
        CgnsBcBoco * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBcRegion );
        int nBcElement = cgnsBcRegion->nElements;
        int nActualBcElement = cgnsBcRegion->GetActualNumberOfBoundaryElements();
        nBFace += nBcElement;
        nActualBcFace += nActualBcElement;

        cout << " iBcRegion  = " << iBcRegion << " numberOfBoundaryElements       = " << nBcElement << "\n";
        cout << " iBcRegion  = " << iBcRegion << " numberOfActualBoundaryElements = " << nActualBcElement << "\n";
    }
    cout << " numberOfBoundaryFaces       = " << nBFace << "\n";
    cout << " numberOfActualBoundaryFaces = " << nActualBcFace << "\n";
    return nActualBcFace;
}

void CgnsZbcBoco::GenerateUnsBcElemConn( CgIntField& bcConn )
{
    int nBcElem = 0;
    int pos = 0;

    cout << " pos = " << pos << "\n";

    for ( int iBcRegion = 0; iBcRegion < this->nBoco; ++ iBcRegion )
    {
        CgnsBcBoco * bcRegion = this->GetCgnsBcRegionBoco( iBcRegion );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );
        SetBcConn( this->cgnsZone, ijkMin, ijkMax, bcConn, pos, nBcElem );
        cout << " pos = " << pos << "\n";
        cout << " nBcElem = " << nBcElem << " boundaryElementSize = " << nBcElem * 4 << "\n";
    }
}


#endif
EndNameSpace