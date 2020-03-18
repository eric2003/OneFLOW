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

#include "CgnsBcConn.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsCoor.h"
#include "CgnsPeriod.h"
#include "CgnsGlobal.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

CgnsBcConn::CgnsBcConn( CgnsZone * cgnsZone )
    : CgnsBcLink( cgnsZone )
{
}

CgnsBcConn::~CgnsBcConn()
{
}

void CgnsBcConn::ReadCgnsBcConnInfo()
{
    int fileId = this->cgnsZone->cgnsBase->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zId = this->cgnsZone->zId;

    CgnsTraits::char33 connName;
    CgnsTraits::char33 donorZoneName;

    cg_conn_info( fileId, baseId, zId, this->bcId,
        connName, & this->gridLocation, & this->gridConnType, & this->pointSetType,
        & nConnPoints, donorZoneName, & donorZoneType, & donorPointSetType, & donorDataType, & nConnDonorPoints );

    this->connName = connName;
    this->donorZoneName  = donorZoneName;

    cout << "\n";
    cout << "   connName      = " << connName << " donorZoneName = " << donorZoneName << "\n";
    cout << "   gridLocation  = " << GridLocationName[ this->gridLocation ] << "\n";
    cout << "   donorDataType = " << DataTypeName[ donorDataType ] << "\n";
    cout << "   gridConnType  = " << GridConnectivityTypeName[ this->gridConnType ] << "\n";
    cout << "   pointSetType  = " << PointSetTypeName[ this->pointSetType ];
    cout << "   donorPointSetType = " << PointSetTypeName[ donorPointSetType ] << "\n";
    cout << "   nConnPoints      = " << nConnPoints << "\n";
    cout << "   nConnDonorPoints = " << nConnDonorPoints << "\n";
}

void CgnsBcConn::ReadCgnsBcConnData()
{
    int fileId = this->cgnsZone->cgnsBase->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zId = this->cgnsZone->zId;

    this->connPoint.resize( nConnPoints );
    this->connDonorPoint.resize( nConnDonorPoints );

    cg_conn_read( fileId, baseId, zId, this->bcId, & this->connPoint[ 0 ], this->donorDataType, & this->connDonorPoint[ 0 ] );
}

void CgnsBcConn::ReadCgnsBcConn()
{
    this->ReadCgnsBcConnInfo();
    this->ReadCgnsBcConnData();
}

void CgnsBcConn::SetPeriodicBc()
{
    CgnsZone * sZone = this->cgnsZone;
    CgnsZone * tZone = ONEFLOW::GetCgnsZoneByName( this->donorZoneName );
    NodeMesh * nodeMesh1 = sZone->cgnsCoor->GetNodeMesh();
    NodeMesh * nodeMesh2 = tZone->cgnsCoor->GetNodeMesh();

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

#endif

EndNameSpace