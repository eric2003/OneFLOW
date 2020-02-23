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

CgnsBcRegionProxy::CgnsBcRegionProxy( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->n1To1     = 0;
    this->nOrdinaryBcRegion = 0;
    this->n1To1General = 0;
    this->nBcRegion = 0;
}

CgnsBcRegionProxy::~CgnsBcRegionProxy()
{
    for ( int ir = 0; ir < this->nOrdinaryBcRegion; ++ ir )
    {
        delete cgnsBcRegions[ ir ];
    }

    for ( int i1To1 = 0; i1To1 < n1To1General; ++ i1To1 )
    {
        delete bcRegion1To1[ i1To1 ];
    }
}


void CgnsBcRegionProxy::CreateCgnsBcRegion()
{
    this->n1To1General = MAX( this->nConn, this->n1To1 );

    cout << "   nConn        = " << this->nConn << endl;
    cout << "   n1To1        = " << this->n1To1 << endl;
    cout << "   n1To1General = " << this->n1To1General << endl;

    this->nBcRegion = this->nOrdinaryBcRegion + this->n1To1General;

    for ( int ir = 0; ir < this->nOrdinaryBcRegion; ++ ir )
    {
        CgnsBcRegion * cgnsBcRegion = new CgnsBcRegion( this->cgnsZone );
        this->AddCgnsBcRegion( cgnsBcRegion );
    }

    for ( int i1To1 = 0; i1To1 < this->n1To1General; ++ i1To1 )
    {
        CgnsBcRegion * cgnsBcRegion = new CgnsBcRegion( this->cgnsZone );
        this->AddCgns1To1BcRegion( cgnsBcRegion );
    }
}

CgnsBcRegion * CgnsBcRegionProxy::GetBcRegion( int ir )
{
    if ( ir < this->nOrdinaryBcRegion )
    {
        return this->cgnsBcRegions[ ir ];
    }
    else
    {
        int irr = ir - this->nOrdinaryBcRegion;
        return this->bcRegion1To1[ irr ];
    }
}

CgnsBcRegion * CgnsBcRegionProxy::GetBcRegion1To1( int i1To1 )
{
    return this->bcRegion1To1[ i1To1 ];
}

void CgnsBcRegionProxy::AddCgnsBcRegion( CgnsBcRegion * cgnsBcRegion )
{
    this->cgnsBcRegions.push_back( cgnsBcRegion );
}

void CgnsBcRegionProxy::AddCgns1To1BcRegion( CgnsBcRegion * cgnsBcRegion )
{
    this->bcRegion1To1.push_back( cgnsBcRegion );
}

void CgnsBcRegionProxy::ConvertToInnerDataStandard()
{
    for ( int ir = 0; ir < nBcRegion; ++ ir )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetBcRegion( ir );
        cgnsBcRegion->ConvertToInnerDataStandard();
    }

    int baseFlag = 1;

    for ( int ir = 0; ir < nBcRegion; ++ ir )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetBcRegion( ir );
        if ( ! cgnsBcRegion->ComputeBase() )
        {
            baseFlag = 0;
            break;
        }
    }

    if ( baseFlag == 0 )
    {
        for ( int ir = 0; ir < nBcRegion; ++ ir )
        {
            CgnsBcRegion * cgnsBcRegion = this->GetBcRegion( ir );
            cgnsBcRegion->ShiftBcRegion();
        }
    }
}

void CgnsBcRegionProxy::ScanBcFace( FaceSolver * face_solver )
{
    cout << " Now ScanBcFace......\n\n";
    cout << " nBcRegion = " << this->nBcRegion << endl;

    for ( int ir = 0; ir < nBcRegion; ++ ir )
    {
        cout << " ir = " << ir << " ";
        CgnsBcRegion * cgnsBcRegion = this->GetBcRegion( ir );
        cout << " BCTypeName = " << ONEFLOW::GetCgnsBcName( cgnsBcRegion->bcType ) << endl;
        cout << " BCRegion Name = " << cgnsBcRegion->name << endl;

        RegionNameMap::AddRegion( cgnsBcRegion->name );
        int bcNameId = RegionNameMap::FindRegionId( cgnsBcRegion->name );
        cgnsBcRegion->nameId = bcNameId;
        cgnsBcRegion->bcId = ir + 1;
        cgnsBcRegion->ScanBcFace( face_solver );
    }
    face_solver->ScanInterfaceBc();
}

void CgnsBcRegionProxy::ReadCgnsGridBoundary()
{
    this->ReadNumberCgnsConnBcInfo();
    this->CreateCgnsBcRegion();

    this->ReadCgnsOrdinaryBcRegion();
    this->ReadCgnsInterfaceBcRegion();

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


void CgnsBcRegionProxy::ReadCgnsInterfaceBcRegion()
{
    int nInterBc = MAX( this->nConn, this->n1To1 );
     for ( int iInterBc = 0; iInterBc < nInterBc; ++ iInterBc )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetBcRegion1To1( iInterBc );

        if ( this->nConn > 0 )
        {
            cgnsBcRegion->ProcessCgns1to1BcRegion( iInterBc + 1 );
        }
        else
        {
            cgnsBcRegion->ReadCgns1to1BoundaryRegion( iInterBc + 1 );
        }
    }
}

void CgnsBcRegionProxy::ReadNumberCgnsConnBcInfo()
{
    this->ReadNumberOfCgnsOrdinaryBcRegions();

    this->ReadNumberOfCgns1To1BcRegions();

    this->ReadNumberOfCgnsConn();
}

void CgnsBcRegionProxy::ReadCgnsOrdinaryBcRegion()
{
    for ( int iBcRegion = 1; iBcRegion <= nOrdinaryBcRegion; ++ iBcRegion )
    {
        cout << "\n-->iBcRegion  = " << iBcRegion;
        cout << " nOrdinaryBcRegion = " << nOrdinaryBcRegion << "\n";
        CgnsBcRegion * cgnsBcRegion = this->GetBcRegion( iBcRegion - 1 );
        cgnsBcRegion->bcId = iBcRegion;
        cgnsBcRegion->ReadCgnsOrdinaryBcRegion();
    }
}


void CgnsBcRegionProxy::ReadNumberOfCgnsOrdinaryBcRegions()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // Determine the number of boundary conditions for this zone.
    cg_nbocos( fileId, baseId, zId, & this->nOrdinaryBcRegion );
}

void CgnsBcRegionProxy::ReadNumberOfCgns1To1BcRegions()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    // find out how many general interfaces there are in this zone
    // the following is the number of structured grid interface
    cg_n1to1( fileId, baseId, zId, & this->n1To1 );

    int n1to1Global = 0;
    cg_n1to1_global( fileId, baseId, & n1to1Global );
}

void CgnsBcRegionProxy::ReadNumberOfCgnsConn()
{
    int fileId = cgnsZone->cgnsBase->fileId;
    int baseId = cgnsZone->cgnsBase->baseId;
    int zId = cgnsZone->zId;

    cg_nconns( fileId, baseId, zId, & this->nConn );
}

void CgnsBcRegionProxy::CreateCgnsBcRegion( CgnsBcRegionProxy * bcRegionProxyIn )
{
    this->nOrdinaryBcRegion = bcRegionProxyIn->nOrdinaryBcRegion;
    this->n1To1             = bcRegionProxyIn->n1To1;
    this->nConn             = bcRegionProxyIn->nConn;

    this->CreateCgnsBcRegion();
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

    for ( int iBcRegion = 0; iBcRegion < this->nBcRegion; ++ iBcRegion )
    {
        CgnsBcRegion * bcRegion = this->GetBcRegion( iBcRegion );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );

        rfact.AddRefBcRegion( ijkMin, ijkMax );
    }

    rfact.Run();

    UInt nnr = rfact.bcregions.size();

    nOrdinaryBcRegion += static_cast<int> (nnr);

    this->nBcRegion = this->nOrdinaryBcRegion + this->n1To1General;

    for ( UInt i = 0; i < nnr; ++ i )
    {
        CgnsBcRegion * rr = new CgnsBcRegion( this->cgnsZone );
        MyRegion * r = rfact.bcregions[ i ];

        int id = static_cast<int> (this->cgnsBcRegions.size() + 1);
        rr->bcId = id;
        rr->ReconstructStrRegion( r->ijkmin, r->ijkmax );

        this->cgnsBcRegions.push_back( rr );
    }
    int kkk = 1;
}

int CgnsBcRegionProxy::GetNumberOfActualBcElements()
{
    int nBFace = 0;
    int nActualBcFace = 0;

    int nBcRegion = this->nBcRegion;

    for ( int iBcRegion = 0; iBcRegion < nBcRegion; ++ iBcRegion )
    {
        CgnsBcRegion * cgnsBcRegion = this->GetBcRegion( iBcRegion );
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

void CgnsBcRegionProxy::GenerateUnsBcElemConn(CgIntField& bcConn )
{
    int nBcElem = 0;
    int pos = 0;

    cout << " pos = " << pos << "\n";

    for ( int iBcRegion = 0; iBcRegion < this->nBcRegion; ++ iBcRegion )
    {
        CgnsBcRegion * bcRegion = this->GetBcRegion( iBcRegion );

        IntField ijkMin( 3 ), ijkMax( 3 );
        bcRegion->ExtractIJKRegionFromBcConn( ijkMin, ijkMax );
        SetBcConn( this->cgnsZone, ijkMin, ijkMax, bcConn, pos, nBcElem );
        cout << " pos = " << pos << "\n";
        cout << " nBcElem = " << nBcElem << " boundaryElementSize = " << nBcElem * 4 << "\n";
    }
}

void CgnsBcRegionProxy::SetPeriodicBc()
{
    for ( int iBcRegion = 0; iBcRegion < this->nBcRegion; ++ iBcRegion )
    {
        CgnsBcRegion * bcRegion = this->GetBcRegion( iBcRegion );
        if ( bcRegion->bcInterface )
        {
            bcRegion->bcInterface->SetPeriodicBc();
        }
    }
}

#endif
EndNameSpace