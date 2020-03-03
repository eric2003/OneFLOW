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

#include "CgnsBcRegionProxy.h"
#include "CgnsBcRegion.h"
#include "CgnsBcInterface.h"
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

CgnsZbcConn::CgnsZbcConn( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->nConn = 0;
}

CgnsZbcConn::~CgnsZbcConn()
{
    for ( int i1To1 = 0; i1To1 < this->nConn; ++ i1To1 )
    {
        delete this->cgnsBcRegionConn[ i1To1 ];
    }
}

void CgnsZbcConn::AddCgnsConnBcRegion( CgnsBcRegion * cgnsBcRegion )
{
    this->cgnsBcRegionConn.push_back( cgnsBcRegion );
}

void CgnsZbcConn::CreateCgnsConnBcRegion()
{
    cout << "   nConn        = " << this->nConn << endl;
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcRegion * cgnsBcRegion = new CgnsBcRegion( this->cgnsZone );
        this->AddCgnsConnBcRegion( cgnsBcRegion );
    }
}

CgnsBcRegionProxy::CgnsBcRegionProxy( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->n1To1 = 0;
    this->nConn = 0;
    this->nBoco = 0;

    this->cgnsZbcConn = new CgnsZbcConn( cgnsZone );
}

CgnsBcRegionProxy::~CgnsBcRegionProxy()
{
    for ( int ir = 0; ir < this->nBoco; ++ ir )
    {
        delete this->cgnsBcRegionBoco[ ir ];
    }

    for ( int i1To1 = 0; i1To1 < this->n1To1; ++ i1To1 )
    {
        delete this->cgnsBcRegion1To1[ i1To1 ];
    }

    for ( int i1To1 = 0; i1To1 < this->nConn; ++ i1To1 )
    {
        delete this->cgnsBcRegionConn[ i1To1 ];
    }

    delete this->cgnsZbcConn;
}

void CgnsBcRegionProxy::CreateCgnsBocoBcRegion()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcRegion * cgnsBcRegion = new CgnsBcRegion( this->cgnsZone );
        this->AddCgnsBocoBcRegion( cgnsBcRegion );
    }
}

void CgnsBcRegionProxy::CreateCgns1To1BcRegion()
{
    cout << "   n1To1        = " << this->n1To1 << endl;

    for ( int i1To1 = 0; i1To1 < this->n1To1; ++ i1To1 )
    {
        CgnsBcRegion * cgnsBcRegion = new CgnsBcRegion( this->cgnsZone );
        this->AddCgns1To1BcRegion( cgnsBcRegion );
    }
}

void CgnsBcRegionProxy::CreateCgnsConnBcRegion()
{
    cout << "   nConn        = " << this->nConn << endl;
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcRegion * cgnsBcRegion = new CgnsBcRegion( this->cgnsZone );
        this->AddCgnsConnBcRegion( cgnsBcRegion );
    }
}

CgnsBcRegion * CgnsBcRegionProxy::GetCgnsBcRegionBoco( int iBoco )
{
    return this->cgnsBcRegionBoco[ iBoco ];
}

CgnsBcRegion * CgnsBcRegionProxy::GetCgnsBcRegion1To1( int i1To1 )
{
    return this->cgnsBcRegion1To1[ i1To1 ];
}

CgnsBcRegion * CgnsBcRegionProxy::GetCgnsBcRegionConn( int iConn )
{
    return this->cgnsBcRegionConn[ iConn ];
}

int CgnsBcRegionProxy::GetNBocoDynamic()
{
    return this->cgnsBcRegionBoco.size();
}

void CgnsBcRegionProxy::AddCgnsBocoBcRegion( CgnsBcRegion * cgnsBcRegion )
{
    this->cgnsBcRegionBoco.push_back( cgnsBcRegion );
}

void CgnsBcRegionProxy::AddCgns1To1BcRegion( CgnsBcRegion * cgnsBcRegion )
{
    this->cgnsBcRegion1To1.push_back( cgnsBcRegion );
}

void CgnsBcRegionProxy::AddCgnsConnBcRegion( CgnsBcRegion * cgnsBcRegion )
{
    this->cgnsBcRegionConn.push_back( cgnsBcRegion );
}

void CgnsBcRegionProxy::ShiftBcRegion()
{
    int baseFlag = 1;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
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
            CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
            cgnsBcRegion->ShiftBcRegion();
        }
    }
}

void CgnsBcRegionProxy::ConvertToInnerDataStandard()
{
    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
        cgnsBcRegion->ConvertToInnerDataStandard();
    }

    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionConn( iConn );
        cgnsBcRegion->ConvertToInnerDataStandard();
    }

    for ( int i1To1 = 0; i1To1 < this->n1To1; ++ i1To1 )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegion1To1( i1To1 );
        cgnsBcRegion->ConvertToInnerDataStandard();
    }

    this->ShiftBcRegion();
}

void CgnsBcRegionProxy::ScanBcFace( FaceSolver * face_solver )
{
    cout << " Now ScanBcFace......\n\n";
    cout << " nBoco = " << this->nBoco << endl;

    for ( int iBoco = 0; iBoco < this->nBoco; ++ iBoco )
    {
        cout << " iBoco = " << iBoco << " ";
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBoco );
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

void CgnsBcRegionProxy::ReadCgnsGridBoundary()
{
    this->ReadCgnsBocoBcRegion();
    this->ReadCgnsConnBcRegion();
    this->ReadCgns1to1BcRegion();

    //this->ReconstructStrRegion();
}

void CgnsBcRegionProxy::FillBcPoints( int * start, int * end, cgsize_t * bcpnts, int dimension )
{
    int icount = 0;
    // lower point of range
    bcpnts[ icount ++ ] = start[ 0 ];
    bcpnts[ icount ++ ] = start[ 1 ];
    if ( dimension == THREE_D )
    {
        bcpnts[ icount ++ ] = start[ 2 ];
    }

    // upper point of range
    bcpnts[ icount ++ ] = end[ 0 ];
    bcpnts[ icount ++ ] = end[ 1 ];
    if ( dimension == THREE_D )
    {
        bcpnts[ icount ++ ] = end[ 2 ];
    }

    cout << " " << start[ 0 ] << " " << end[ 0 ];
    cout << " " << start[ 1 ] << " " << end[ 1 ];
    if ( dimension == THREE_D )
    {
        cout << " " << start[ 2 ] << " " << end[ 2 ];
    }
    cout << "\n";
}

void CgnsBcRegionProxy::FillBcPoints3D( int * start, int * end, cgsize_t * bcpnts )
{
    int icount = 0;
    // lower point of range
    bcpnts[ icount ++ ] = start[ 0 ];
    bcpnts[ icount ++ ] = start[ 1 ];
    bcpnts[ icount ++ ] = start[ 2 ];

    // upper point of range
    bcpnts[ icount ++ ] = end[ 0 ];
    bcpnts[ icount ++ ] = end[ 1 ];
    bcpnts[ icount ++ ] = end[ 2 ];

    cout << " " << start[ 0 ] << " " << end[ 0 ];
    cout << " " << start[ 1 ] << " " << end[ 1 ];
    cout << " " << start[ 2 ] << " " << end[ 2 ];
    cout << "\n";
}

void CgnsBcRegionProxy::FillRegion( TestRegion * r, cgsize_t * ipnts, int dimension )
{
    //int dimension = cgnsZone->cgnsBase->celldim;
    int icount = 0;
    //lower point of receiver range
    ipnts[ icount ++ ] = r->p1[ 0 ];
    ipnts[ icount ++ ] = r->p1[ 1 ];
    if ( dimension == THREE_D )
    {
        ipnts[ icount ++ ] = r->p1[ 2 ];
    }
    //upper point of receiver range
    ipnts[ icount ++ ] = r->p2[ 0 ];
    ipnts[ icount ++ ] = r->p2[ 1 ];
    if ( dimension == THREE_D )
    {
        ipnts[ icount ++ ] = r->p2[ 2 ];
    }
}

void CgnsBcRegionProxy::FillInterface( BcRegion * bcRegion, cgsize_t * ipnts, cgsize_t * ipntsdonor, int * itranfrm, int dimension )
{
    TestRegionM trm;
    trm.Run( bcRegion, dimension );
    this->FillRegion( & trm.s, ipnts, dimension );
    this->FillRegion( & trm.t, ipntsdonor, dimension );

    // set up Transform
    itranfrm[ 0 ] = trm.itransform[ 0 ];
    itranfrm[ 1 ] = trm.itransform[ 1 ];
    itranfrm[ 2 ] = trm.itransform[ 2 ];
    cout << " itranfrm = ";
    cout << itranfrm[ 0 ] << " ";
    cout << itranfrm[ 1 ] << " ";
    cout << itranfrm[ 2 ] << " ";
    cout << "\n";
}

void CgnsBcRegionProxy::DumpCgnsGridBoundary( Grid * gridIn )
{
    StrGrid * grid = StrGridCast( gridIn );

    BcRegionGroup * bcRegionGroup = grid->bcRegionGroup;

    int nBcRegions = bcRegionGroup->regions->size();

    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zoneId = cgnsZone->zId;

    cout << " fildId = " << fileId << " baseId = " << baseId << " zoneId = " << zoneId << "\n";

    BcTypeMap * bcTypeMap = new BcTypeMap();
    bcTypeMap->Init();

    cgsize_t ipnts[ 6 ], ipntsdonor[ 6 ];
    int itranfrm[ 3 ];

    for ( int ir = 0; ir < nBcRegions; ++ ir )
    {
        BcRegion * bcRegion = bcRegionGroup->GetBcRegion( ir );

        BCType_t bctype = static_cast< BCType_t >( bcTypeMap->OneFlow2Cgns( bcRegion->bcType ) );
        int dimension = cgnsZone->cgnsBase->celldim;
        if ( bctype == BCTypeNull )
        {
            FillInterface( bcRegion, ipnts, ipntsdonor, itranfrm, dimension );
            int zid = bcRegion->t->zid - 1;
            Grid * tGrid = GlobalGrid::GetGrid( zid );
            string & donorName = tGrid->name;
            // write 1-to-1 info
            int index_conn = -1;
            cg_1to1_write( fileId, baseId, zoneId, bcRegion->regionName.c_str(), donorName.c_str(), ipnts, ipntsdonor,itranfrm, & index_conn );
            cout << " regionName = " << bcRegion->regionName << " donorName = " << donorName << " index_conn = " << index_conn << "\n";
        }
        else
        {
            BasicRegion * s = bcRegion->s;
            FillBcPoints( s->start, s->end, ipnts, dimension );
            //FillBcPoints3D( s->start, s->end, ipnts );
            int bcId = -1;
            cg_boco_write( fileId, baseId, zoneId, bcRegion->regionName.c_str(), bctype, PointRange, 2, ipnts, &bcId );
            cout << " bcId = " << bcId << " regionName = " << bcRegion->regionName << "\n";
        }
    }

    delete bcTypeMap;
}

void CgnsBcRegionProxy::ReadNumberOfCgnsConn()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    cg_nconns( fileId, baseId, zId, & this->nConn );
}

void CgnsBcRegionProxy::ReadNumberOfCgns1To1()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // find out how many general interfaces there are in this zone
    // the following is the number of structured grid interface
    cg_n1to1( fileId, baseId, zId, & this->n1To1 );
}

void CgnsBcRegionProxy::ReadCgns1to1BcRegion()
{
    this->ReadNumberOfCgns1To1();
    this->CreateCgns1To1BcRegion();

    for ( int i1To1 = 0; i1To1 < this->n1To1; ++ i1To1 )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegion1To1( i1To1 );
        cgnsBcRegion->ReadCgns1to1BcRegion( i1To1 + 1 );
    }
}

void CgnsBcRegionProxy::ReadCgnsConnBcRegion()
{
    this->ReadNumberOfCgnsConn();
    this->CreateCgnsConnBcRegion();
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionConn( iConn );
        cgnsBcRegion->ReadCgnsConnBcRegion( iConn + 1 );
    }
}

void CgnsBcRegionProxy::ReadCgnsBocoBcRegion()
{
    this->ReadNumberOfCgnsBoco();
    this->CreateCgnsBocoBcRegion();

    for ( int iBcRegion = 0; iBcRegion < nBoco; ++ iBcRegion )
    {
        cout << "\n-->iBcRegion  = " << iBcRegion;
        cout << " nOrdinaryBcRegion = " << nBoco << "\n";
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBcRegion );
        cgnsBcRegion->bcId = iBcRegion + 1;
        cgnsBcRegion->ReadCgnsBocoBcRegion();
    }
}

void CgnsBcRegionProxy::ReadNumberOfCgnsBoco()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // Determine the number of boundary conditions for this zone.
    cg_nbocos( fileId, baseId, zId, & this->nBoco );
}

void CgnsBcRegionProxy::CreateCgnsBcRegion( CgnsBcRegionProxy * bcRegionProxyIn )
{
    this->nBoco = bcRegionProxyIn->nBoco;
    this->CreateCgnsBocoBcRegion();

    this->n1To1 = bcRegionProxyIn->n1To1;
    this->CreateCgns1To1BcRegion();

    this->nConn = bcRegionProxyIn->nConn;
    this->CreateCgnsConnBcRegion();
}

void CgnsBcRegionProxy::ReconstructStrRegion()
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
        CgnsBcRegion * bcRegion = this->GetCgnsBcRegionBoco( iBoco );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );

        rfact.AddRefBcRegion( ijkMin, ijkMax );
    }

    rfact.Run();

    UInt nnr = rfact.bcregions.size();

    nBoco += static_cast<int> (nnr);

    for ( UInt i = 0; i < nnr; ++ i )
    {
        CgnsBcRegion * rr = new CgnsBcRegion( this->cgnsZone );
        MyRegion * r = rfact.bcregions[ i ];

        int id = static_cast<int> ( this->GetNBocoDynamic() + 1 );
        rr->bcId = id;
        rr->ReconstructStrRegion( r->ijkmin, r->ijkmax );

        this->AddCgnsBocoBcRegion( rr );
    }
    int kkk = 1;
}

int CgnsBcRegionProxy::GetNumberOfActualBcElements()
{
    int nBFace = 0;
    int nActualBcFace = 0;

    for ( int iBcRegion = 0; iBcRegion < this->nBoco; ++ iBcRegion )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionBoco( iBcRegion );
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

void CgnsBcRegionProxy::GenerateUnsBcElemConn( CgIntField& bcConn )
{
    int nBcElem = 0;
    int pos = 0;

    cout << " pos = " << pos << "\n";

    for ( int iBcRegion = 0; iBcRegion < this->nBoco; ++ iBcRegion )
    {
        CgnsBcRegion * bcRegion = this->GetCgnsBcRegionBoco( iBcRegion );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );
        SetBcConn( this->cgnsZone, ijkMin, ijkMax, bcConn, pos, nBcElem );
        cout << " pos = " << pos << "\n";
        cout << " nBcElem = " << nBcElem << " boundaryElementSize = " << nBcElem * 4 << "\n";
    }
}

void CgnsBcRegionProxy::SetPeriodicBc()
{
    for ( int iConn = 0; iConn < this->nConn; ++ iConn )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegionConn( iConn );
        cgnsBcRegion->bcInterface->SetPeriodicBc();
    }

    for ( int i1To1 = 0; i1To1 < this->n1To1; ++ i1To1 )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetCgnsBcRegion1To1( i1To1 );
        cgnsBcRegion->bcInterface->SetPeriodicBc();
    }
}

#endif
EndNameSpace