#include "GridHandles.h"
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

#include "CgnsZbc.h"
#include "CgnsBcBoco.h"
#include "CgnsZbc1to1.h"
#include "CgnsZbcConn.h"
#include "CgnsZbcBoco.h"
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

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsZbc::CgnsZbc( CgnsZone & cgnsZone )
    : cgnsZbcConn( std::make_unique< CgnsZbcConn >( cgnsZone ) ),
      cgnsZbc1to1( std::make_unique< CgnsZbc1to1 >( cgnsZone ) ),
      cgnsZbcBoco( std::make_unique< CgnsZbcBoco >( cgnsZone ) ),
      cgnsZone( cgnsZone )
{
}

CgnsZbc::~CgnsZbc() = default;

CgnsZbcConn & CgnsZbc::RequireCgnsZbcConn()
{
    if ( this->cgnsZbcConn == nullptr )
    {
        throw std::logic_error( "CgnsZbc: CgnsZbcConn is not initialized" );
    }
    return *this->cgnsZbcConn;
}

CgnsZbc1to1 & CgnsZbc::RequireCgnsZbc1to1()
{
    if ( this->cgnsZbc1to1 == nullptr )
    {
        throw std::logic_error( "CgnsZbc: CgnsZbc1to1 is not initialized" );
    }
    return *this->cgnsZbc1to1;
}

CgnsZbcBoco & CgnsZbc::RequireCgnsZbcBoco()
{
    if ( this->cgnsZbcBoco == nullptr )
    {
        throw std::logic_error( "CgnsZbc: CgnsZbcBoco is not initialized" );
    }
    return *this->cgnsZbcBoco;
}

void CgnsZbc::ConvertToInnerDataStandard()
{
    this->RequireCgnsZbcBoco().ConvertToInnerDataStandard();

    this->RequireCgnsZbcConn().ConvertToInnerDataStandard();

    this->RequireCgnsZbc1to1().ConvertToInnerDataStandard();

    this->RequireCgnsZbcBoco().ShiftBcRegion();
}

void CgnsZbc::ScanBcFace( FaceSolver & faceSolver )
{
    this->RequireCgnsZbcBoco().ScanBcFace( faceSolver );
}

void CgnsZbc::ReadCgnsGridBoundary()
{
    this->RequireCgnsZbcBoco().ReadCgnsZbcBoco();
    this->RequireCgnsZbcConn().ReadCgnsZbcConn();
    this->RequireCgnsZbc1to1().ReadCgnsZbc1to1();
}

void CgnsZbc::DumpCgnsGridBoundary()
{
    this->RequireCgnsZbcBoco().DumpCgnsZbcBoco();
    this->RequireCgnsZbcConn().DumpCgnsZbcConn();
    this->RequireCgnsZbc1to1().DumpCgnsZbc1to1();
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

    std::cout << " " << start[ 0 ] << " " << end[ 0 ];
    std::cout << " " << start[ 1 ] << " " << end[ 1 ];
    if ( dimension == THREE_D )
    {
        std::cout << " " << start[ 2 ] << " " << end[ 2 ];
    }
    std::cout << "\n";
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

    std::cout << " " << start[ 0 ] << " " << end[ 0 ];
    std::cout << " " << start[ 1 ] << " " << end[ 1 ];
    std::cout << " " << start[ 2 ] << " " << end[ 2 ];
    std::cout << "\n";
}

void CgnsZbc::FillRegion( TestRegion * r, cgsize_t * ipnts, int dimension )
{
    //int dimension = cgnsZone.cgnsBase.celldim;
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

    // std::set up Transform
    itranfrm[ 0 ] = trm.itransform[ 0 ];
    itranfrm[ 1 ] = trm.itransform[ 1 ];
    itranfrm[ 2 ] = trm.itransform[ 2 ];
    std::cout << " itranfrm = ";
    std::cout << itranfrm[ 0 ] << " ";
    std::cout << itranfrm[ 1 ] << " ";
    std::cout << itranfrm[ 2 ] << " ";
    std::cout << "\n";
}

void CgnsZbc::DumpCgnsGridBoundary( Grid * gridIn, const Grids & grids )
{
    StrGrid * grid = StrGridCast( gridIn );

    BcRegionGroup * bcRegionGroup = grid->bcRegionGroup.get();

    int nBcRegions = bcRegionGroup->regions.size();

    int fileId = cgnsZone.cgnsBase.cgnsFile->fileId;
    int baseId = cgnsZone.cgnsBase.baseId;
    int zoneId = cgnsZone.zId;

    std::cout << " fildId = " << fileId << " baseId = " << baseId << " zoneId = " << zoneId << "\n";

    BcTypeMap bcTypeMap;
    bcTypeMap.Init();

    cgsize_t ipnts[ 6 ], ipntsdonor[ 6 ];
    int itranfrm[ 3 ];

    for ( int ir = 0; ir < nBcRegions; ++ ir )
    {
        BcRegion * bcRegion = bcRegionGroup->GetBcRegion( ir );

        BCType_t bctype = static_cast< BCType_t >( bcTypeMap.OneFlow2Cgns( bcRegion->bcType ) );
        int dimension = cgnsZone.cgnsBase.celldim;
        if ( bctype == BCTypeNull )
        {
            FillInterface( bcRegion, ipnts, ipntsdonor, itranfrm, dimension );
            int zid = bcRegion->t->zid - 1;
            const Grid & tGrid = GridAt( grids, zid );
            const std::string & donorName = tGrid.name;
            // write 1-to-1 info
            int index_conn = -1;
            cg_1to1_write( fileId, baseId, zoneId, bcRegion->regionName.c_str(), donorName.c_str(), ipnts, ipntsdonor,itranfrm, & index_conn );
            std::cout << " regionName = " << bcRegion->regionName << " donorName = " << donorName << " index_conn = " << index_conn << "\n";
        }
        else
        {
            BasicRegion * s = bcRegion->s.get();
            FillBcPoints( s->start, s->end, ipnts, dimension );
            //FillBcPoints3D( s->start, s->end, ipnts );
            int bcId = -1;
            cg_boco_write( fileId, baseId, zoneId, bcRegion->regionName.c_str(), bctype, PointRange, 2, ipnts, &bcId );
            std::cout << " bcId = " << bcId << " regionName = " << bcRegion->regionName << "\n";
        }
    }

}

void CgnsZbc::CreateCgnsZbc( CgnsZbc * cgnsZbcIn )
{
    this->RequireCgnsZbcBoco().ReadZnboco( cgnsZbcIn->RequireCgnsZbcBoco().nBoco );
    this->RequireCgnsZbcBoco().CreateCgnsZbc();

    this->RequireCgnsZbc1to1().ReadZn1to1( cgnsZbcIn->RequireCgnsZbc1to1().n1to1 );
    this->RequireCgnsZbc1to1().CreateCgnsZbc();

    this->RequireCgnsZbcConn().ReadZnconn( cgnsZbcIn->RequireCgnsZbcConn().nConn );
    this->RequireCgnsZbcConn().CreateCgnsZbc();
}

int CgnsZbc::GetNumberOfActualBcElements()
{
    return this->RequireCgnsZbcBoco().GetNumberOfActualBcElements();
}

void CgnsZbc::GenerateUnsBcElemConn( CgIntField& bcConn )
{
    this->RequireCgnsZbcBoco().GenerateUnsBcElemConn( bcConn );
}

void CgnsZbc::SetPeriodicBc()
{
    this->RequireCgnsZbcConn().SetPeriodicBc();

    this->RequireCgnsZbc1to1().SetPeriodicBc();
}

#endif
EndNameSpace
