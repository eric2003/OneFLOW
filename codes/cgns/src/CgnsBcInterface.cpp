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

#include "CgnsBcInterface.h"
#include "CgnsBcRegion.h"
#include "CgnsBc1to1.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsPeriod.h"
#include "CgnsGlobal.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS


CgnsBcInterface::CgnsBcInterface( CgnsBcRegion * bcRegion )
{
    this->bcRegion = bcRegion;
    this->flag1To1 = false;
}

CgnsBcInterface::~CgnsBcInterface()
{
    ;
}

void CgnsBcInterface::ReadCgnsBcConnInfo()
{
    int fileId = this->bcRegion->cgnsZone->cgnsBase->fileId;
    int baseId = this->bcRegion->cgnsZone->cgnsBase->baseId;
    int zId = this->bcRegion->cgnsZone->zId;

    CgnsTraits::char33 connName;
    CgnsTraits::char33 donorZoneName;

    cg_conn_info( fileId, baseId, zId, this->bcRegion->bcId,
                  connName, & this->bcRegion->gridLocation, & this->bcRegion->gridConnType, & this->bcRegion->pointSetType,
                  & nConnPoints, donorZoneName, & donorZoneType, & donorPointSetType, & donorDataType, & nConnDonorPoints );

    this->bcRegion->name = connName;
    this->donorZoneName  = donorZoneName;

    cout << "\n";
    cout << "   connName      = " << connName << " donorZoneName = " << donorZoneName << "\n";
    cout << "   gridLocation  = " << GridLocationName[ this->bcRegion->gridLocation ] << "\n";
    cout << "   donorDataType = " << DataTypeName[ donorDataType ] << "\n";
    cout << "   gridConnType  = " << GridConnectivityTypeName[ this->bcRegion->gridConnType ] << "\n";
    cout << "   pointSetType  = " << PointSetTypeName[ this->bcRegion->pointSetType ];
    cout << "   donorPointSetType = " << PointSetTypeName[ donorPointSetType ] << "\n";
    cout << "   nConnPoints      = " << nConnPoints << "\n";
    cout << "   nConnDonorPoints = " << nConnDonorPoints << "\n";
}

void CgnsBcInterface::ReadCgnsBcConnData()
{
    int fileId = this->bcRegion->cgnsZone->cgnsBase->fileId;
    int baseId = this->bcRegion->cgnsZone->cgnsBase->baseId;
    int zId = this->bcRegion->cgnsZone->zId;

    this->connPoint.resize( nConnPoints );
    this->connDonorPoint.resize( nConnDonorPoints );

    cg_conn_read( fileId, baseId, zId, this->bcRegion->bcId, & this->connPoint[ 0 ], this->donorDataType, & this->connDonorPoint[ 0 ] );
}


void CgnsBcInterface::SetPeriodicBc()
{
    if ( this->flag1To1 ) return;

    CgnsZone * sZone = this->bcRegion->cgnsZone;
    CgnsZone * tZone = ONEFLOW::GetCgnsZoneByName( this->donorZoneName );
    NodeMesh * nodeMesh1 = sZone->nodeMesh;
    NodeMesh * nodeMesh2 = tZone->nodeMesh;

    for ( int i = 0; i < nConnPoints; ++ i )
    {
        if ( i == 40 && tZone->zId == 30  )
        {
            int kkk = 1;
        }
        int id1 = this->connPoint[ i ];
        int id2 = this->connDonorPoint[ i ];

        CgIntField fNodeId1, fNodeId2;
        sZone->GetElementNodeId( id1, fNodeId1 );
        tZone->GetElementNodeId( id2, fNodeId2 );

        f2fmap.AddFacePoint(fNodeId1, fNodeId2, nodeMesh1, nodeMesh2 );
        int kkk = 1;
    }
}

void CgnsBcInterface::ReadCgnsBc1To1()
{
    int fileId = this->bcRegion->cgnsZone->cgnsBase->fileId;
    int baseId = this->bcRegion->cgnsZone->cgnsBase->baseId;
    int zId = this->bcRegion->cgnsZone->zId;

    this->flag1To1 = true;

    this->nConnPoints = 6;
    this->nConnDonorPoints = 6;

    this->connPoint.resize( nConnPoints );
    this->connDonorPoint.resize( nConnDonorPoints );

    this->bcRegion->gridConnType = CGNS_ENUMV( Abutting1to1 );
    this->bcRegion->bcType       = CGNS_ENUMV( BCTypeNull );
    this->bcRegion->pointSetType = CGNS_ENUMV( PointRange );
    this->bcRegion->gridLocation = CGNS_ENUMV( FaceCenter );

    this->donorPointSetType = CGNS_ENUMV( PointRange );
    this->donorDataType     = CGNS_ENUMV( Integer );

    CgnsTraits::char33 connName;
    CgnsTraits::char33 donorZoneName;

    //Zone Connectivity
    cg_goto( fileId, baseId, "Zone_t", zId, "ZoneGridConnectivity_t", 1, "GridConnectivity1to1_t", 1, "end" );

    cg_1to1_read( fileId, baseId, zId, this->bcRegion->bcId, connName, donorZoneName, & this->connPoint[ 0 ], & this->connDonorPoint[ 0 ], itranfrm );

    this->bcRegion->name = connName;
    this->donorZoneName  = donorZoneName;

    cout << "\n";
    cout << "   connName      = " << connName << " donorZoneName = " << donorZoneName << "\n";
    cout << "   gridLocation  = " << GridLocationName[ this->bcRegion->gridLocation ] << "\n";
    cout << "   donorDataType = " << DataTypeName[ this->donorDataType ] << "\n";
    cout << "   gridConnType  = " << GridConnectivityTypeName[ this->bcRegion->gridConnType ] << "\n";
    cout << "   pointSetType  = " << PointSetTypeName[ this->bcRegion->pointSetType ];
    cout << "   donorPointSetType = " << PointSetTypeName[ this->donorPointSetType ] << "\n";
    cout << "   nConnPoints      = " << nConnPoints << "\n";
    cout << "   nConnDonorPoints = " << nConnDonorPoints << "\n";

    cout << "   range (this zone )= ";
    int width = 5;
    cout << setw( width ) << connPoint[ 0 ];
    cout << setw( width ) << connPoint[ 1 ];
    cout << setw( width ) << connPoint[ 2 ] << "\n";
    cout << "                       ";
    cout << setw( width ) << connPoint[ 3 ];
    cout << setw( width ) << connPoint[ 4 ];
    cout << setw( width ) << connPoint[ 5 ] << "\n";
    cout << "   range (donor zone)= ";
    cout << setw( width ) << connDonorPoint[ 0 ];
    cout << setw( width ) << connDonorPoint[ 1 ];
    cout << setw( width ) << connDonorPoint[ 2 ] << "\n";
    cout << "                       ";
    cout << setw( width ) << connDonorPoint[ 3 ];
    cout << setw( width ) << connDonorPoint[ 4 ];
    cout << setw( width ) << connDonorPoint[ 5 ] << "\n";
    cout << "   transform = " << itranfrm[ 0 ] << " " << itranfrm[ 1 ] << " " << itranfrm[ 2 ] << "\n";

    int transform[ 3 ][ 3 ];

    //For 3-D, the
    //transformation matrix T is constructed from Transform = [¡Àa, ¡Àb, ¡Àc ] as follows:
    //T =
    //[ sgn( a ) del( a ? 1 ) sgn( b ) del( b ? 1 ) sgn( c ) del( c ? 1 )]
    //[ sgn( a ) del( a ? 2 ) sgn( b ) del( b ? 2 ) sgn( c ) del( c ? 2 )]
    //[ sgn( a ) del( a ? 3 ) sgn( b ) del( b ? 3 ) sgn( c ) del( c ? 3 )]

    int celldim = this->bcRegion->cgnsZone->cgnsBase->celldim;

    for ( int i = 0; i < celldim; ++ i )
    {
        for ( int j = 0; j < celldim; ++ j )
        {
            transform[ i ][ j ] = SIGN( 1, itranfrm[ j ] ) * AbsoluteDiagonalId( itranfrm[ j ], i + 1 );
        }
    }
}

void CgnsBcInterface::ConvertToInnerDataStandard()
{
    for ( int eId = 0; eId < this->nConnPoints; ++ eId )
    {
        this->connPoint[ eId ] -= 1;
    }

    for ( int eId = 0; eId < this->nConnDonorPoints; ++ eId )
    {
        this->connDonorPoint[ eId ] -= 1;
    }
}

void CgnsBcInterface::ShiftBcRegion()
{
    int nTCell = this->bcRegion->cgnsZone->nCell;
    if ( this->bcRegion->gridLocation != Vertex )
    {
        for ( int eId = 0; eId < this->nConnPoints; ++ eId )
        {
            this->connPoint[ eId ] += nTCell;
        }

        for ( int eId = 0; eId < this->nConnDonorPoints; ++ eId )
        {
            this->connDonorPoint[ eId ] += nTCell;
        }
    }
}
#endif
EndNameSpace