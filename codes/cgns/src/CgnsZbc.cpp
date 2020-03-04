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

#include "CgnsZbc.h"
#include "CgnsBcBoco.h"
#include "CgnsZbc1to1.h"
#include "CgnsZbcConn.h"
#include "CgnsZbcBoco.h"
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

CgnsZbc::CgnsZbc( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;

    this->cgnsZbcConn = new CgnsZbcConn( cgnsZone );
    this->cgnsZbc1to1 = new CgnsZbc1to1( cgnsZone );
    this->cgnsZbcBoco = new CgnsZbcBoco( cgnsZone );

}

CgnsZbc::~CgnsZbc()
{
    delete this->cgnsZbcConn;
    delete this->cgnsZbc1to1;
    delete this->cgnsZbcBoco;
}

void CgnsZbc::ConvertToInnerDataStandard()
{
    this->cgnsZbcBoco->ConvertToInnerDataStandard();

    this->cgnsZbcConn->ConvertToInnerDataStandard();

    this->cgnsZbc1to1->ConvertToInnerDataStandard();

    this->cgnsZbcBoco->ShiftBcRegion();
}

void CgnsZbc::ScanBcFace( FaceSolver * face_solver )
{
    this->cgnsZbcBoco->ScanBcFace( face_solver );
}

void CgnsZbc::ReadCgnsGridBoundary()
{
    this->cgnsZbcBoco->ReadCgnsZbcBoco();
    this->cgnsZbcConn->ReadCgnsZbcConn();
    this->cgnsZbc1to1->ReadCgnsZbc1to1();
}

void CgnsZbc::FillBcPoints( int * start, int * end, cgsize_t * bcpnts, int dimension )
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

void CgnsZbc::FillBcPoints3D( int * start, int * end, cgsize_t * bcpnts )
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

void CgnsZbc::FillRegion( TestRegion * r, cgsize_t * ipnts, int dimension )
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

void CgnsZbc::FillInterface( BcRegion * bcRegion, cgsize_t * ipnts, cgsize_t * ipntsdonor, int * itranfrm, int dimension )
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

void CgnsZbc::DumpCgnsGridBoundary( Grid * gridIn )
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

void CgnsZbc::CreateCgnsZbc( CgnsZbc * cgnsZbcIn )
{
    this->cgnsZbcBoco->ReadZnboco( cgnsZbcIn->cgnsZbcBoco->nBoco );
    this->cgnsZbcBoco->CreateCgnsZbc();

    this->cgnsZbc1to1->ReadZn1to1( cgnsZbcIn->cgnsZbc1to1->n1to1 );
    this->cgnsZbc1to1->CreateCgnsZbc();

    this->cgnsZbcConn->ReadZnconn( cgnsZbcIn->cgnsZbcConn->nConn );
    this->cgnsZbcConn->CreateCgnsZbc();
}

int CgnsZbc::GetNumberOfActualBcElements()
{
    return this->cgnsZbcBoco->GetNumberOfActualBcElements();
}

void CgnsZbc::GenerateUnsBcElemConn( CgIntField& bcConn )
{
    this->cgnsZbcBoco->GenerateUnsBcElemConn( bcConn );
}

void CgnsZbc::SetPeriodicBc()
{
    this->cgnsZbcConn->SetPeriodicBc();

    this->cgnsZbc1to1->SetPeriodicBc();
}

#endif
EndNameSpace