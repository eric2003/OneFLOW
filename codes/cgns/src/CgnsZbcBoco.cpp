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
        delete this->cgnsBcBocos[ iBoco ];
    }
}

void CgnsZbcBoco::AddCgnsBcBoco( CgnsBcBoco * cgnsBcBoco )
{
    this->cgnsBcBocos.push_back( cgnsBcBoco );
    int id = this->cgnsBcBocos.size();
    cgnsBcBoco->bcId = id;
}

CgnsBcBoco * CgnsZbcBoco::GetCgnsBc( int iBoco )
{
    return this->cgnsBcBocos[ iBoco ];
}

void CgnsZbcBoco::CreateCgnsZbc()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcBoco = new CgnsBcBoco( this->cgnsZone );
        this->AddCgnsBcBoco( cgnsBcBoco );
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

void CgnsZbcBoco::ScanBcFace( FaceSolver * face_solver )
{
    cout << " Now ScanBcFace......\n\n";
    cout << " nBoco = " << this->nBoco << endl;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        cout << " iBoco = " << iBoco << " ";
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        cout << " BCTypeName = " << ONEFLOW::GetCgnsBcName( cgnsBcBoco->bcType ) << endl;
        cout << " BCRegion Name = " << cgnsBcBoco->name << endl;

        RegionNameMap::AddRegion( cgnsBcBoco->name );
        int bcNameId = RegionNameMap::FindRegionId( cgnsBcBoco->name );
        cgnsBcBoco->nameId = bcNameId;
        cgnsBcBoco->ScanBcFace( face_solver );
    }
    face_solver->ScanInterfaceBc();
}

void CgnsZbcBoco::PrintZnboco()
{
    cout << "   nBoco        = " << this->nBoco << endl;
}

void CgnsZbcBoco::ReadZnboco()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

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
        cout << "\n";
        cout << "-->iBoco  = " << iBoco << " nBoco = " << nBoco << "\n";
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        cgnsBcBoco->ReadCgnsBcBoco();
    }
}

int CgnsZbcBoco::GetNumberOfActualBcElements()
{
    int nBFace = 0;
    int nActualBcFace = 0;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * cgnsBcBoco = this->GetCgnsBc( iBoco );
        int nBcElement = cgnsBcBoco->nElements;
        int nActualBcElement = cgnsBcBoco->GetActualNumberOfBoundaryElements();
        nBFace += nBcElement;
        nActualBcFace += nActualBcElement;

        cout << " iBoco  = " << iBoco << " numberOfBoundaryElements       = " << nBcElement << "\n";
        cout << " iBoco  = " << iBoco << " numberOfActualBoundaryElements = " << nActualBcElement << "\n";
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

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcBoco * bcRegion = this->GetCgnsBc( iBoco );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );
        SetBcConn( this->cgnsZone, ijkMin, ijkMax, bcConn, pos, nBcElem );
        cout << " pos = " << pos << "\n";
        cout << " nBcElem = " << nBcElem << " boundaryElementSize = " << nBcElem * 4 << "\n";
    }
}


#endif
EndNameSpace