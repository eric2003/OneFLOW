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

#include "CgnsBcConn.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "CgnsCoor.h"
#include "CgnsFile.h"
#include "CgnsPeriod.h"
#include "CgnsGlobal.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>
#include <iomanip>

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
    int fileId = this->cgnsZone->cgnsBase.cgnsFile->fileId;
    int baseId = this->cgnsZone->cgnsBase.baseId;
    int zId = this->cgnsZone->zId;

    CgnsTraits::char33 connName;
    CgnsTraits::char33 donorZoneName;

    CgInt nConnPoints = 0;
    CgInt nConnDonorPoints = 0;

    cg_conn_info( fileId, baseId, zId, this->bcId,
        connName, & this->gridLocation, & this->gridConnType, & this->pointSetType,
        & nConnPoints, donorZoneName, & donorZoneType, & donorPointSetType, & donorDataType, & nConnDonorPoints );

    this->SetNConnPoints( nConnPoints );
    this->SetNConnDonorPoints( nConnDonorPoints );

    this->connName = connName;
    this->donorZoneName  = donorZoneName;

    std::cout << "\n";
    std::cout << "   connName      = " << connName << " donorZoneName = " << donorZoneName << "\n";
    std::cout << "   gridLocation  = " << GridLocationName[ this->gridLocation ] << "\n";
    std::cout << "   donorDataType = " << DataTypeName[ donorDataType ] << "\n";
    std::cout << "   gridConnType  = " << GridConnectivityTypeName[ this->gridConnType ] << "\n";
    std::cout << "   pointSetType  = " << PointSetTypeName[ this->pointSetType ];
    std::cout << "   donorPointSetType = " << PointSetTypeName[ donorPointSetType ] << "\n";
    std::cout << "   nConnPoints      = " << this->GetNConnPoints() << "\n";
    std::cout << "   nConnDonorPoints = " << this->GetNConnDonorPoints() << "\n";
}

void CgnsBcConn::DumpCgnsBcConnInfo()
{
    std::cout << "\n";
    std::cout << "   connName      = " << this->connName
        << " donorZoneName = " << this->donorZoneName << "\n";
    std::cout << "   gridLocation  = "
        << GridLocationName[ this->gridLocation ] << "\n";
    std::cout << "   donorDataType = "
        << DataTypeName[ this->donorDataType ] << "\n";
    std::cout << "   gridConnType  = "
        << GridConnectivityTypeName[ this->gridConnType ] << "\n";
    std::cout << "   pointSetType  = "
        << PointSetTypeName[ this->pointSetType ];
    std::cout << "   donorPointSetType = "
        << PointSetTypeName[ this->donorPointSetType ] << "\n";
    std::cout << "   nConnPoints      = "
        << this->GetNConnPoints() << "\n";
    std::cout << "   nConnDonorPoints = "
        << this->GetNConnDonorPoints() << "\n";
}

void CgnsBcConn::ReadCgnsBcConnData()
{
    int fileId = this->cgnsZone->cgnsBase.cgnsFile->fileId;
    int baseId = this->cgnsZone->cgnsBase.baseId;
    int zId = this->cgnsZone->zId;

    this->ResizeConnPoint( this->GetNConnPoints() );
    this->ResizeConnDonorPoint( this->GetNConnDonorPoints() );

    cg_conn_read( fileId, baseId, zId, this->bcId, this->GetConnPointData(), this->donorDataType, this->GetConnDonorPointData() );
}

void CgnsBcConn::DumpCgnsBcConnData()
{
    int fileId = this->cgnsZone->cgnsBase.cgnsFile->fileId;
    int baseId = this->cgnsZone->cgnsBase.baseId;
    int zId = this->cgnsZone->zId;

    //this->connPoint.resize( nConnPoints );
    //this->connDonorPoint.resize( nConnDonorPoints );

    //cg_conn_read( fileId, baseId, zId, this->bcId, & this->connPoint[ 0 ], this->donorDataType, & this->connDonorPoint[ 0 ] );
}

void CgnsBcConn::ReadCgnsBcConn()
{
    this->ReadCgnsBcConnInfo();
    this->ReadCgnsBcConnData();
}

void CgnsBcConn::DumpCgnsBcConn()
{
    this->DumpCgnsBcConnInfo();
    this->DumpCgnsBcConnData();
}

void CgnsBcConn::SetPeriodicBc()
{
    CgnsZone * sZone = this->cgnsZone;
    CgnsZone * tZone = ONEFLOW::GetCgnsZoneByName( this->donorZoneName );
    NodeMesh * nodeMesh1 = sZone->RequireCgnsCoor().GetNodeMesh();
    NodeMesh * nodeMesh2 = tZone->RequireCgnsCoor().GetNodeMesh();

    for ( int i = 0; i < this->GetNConnPoints(); ++ i )
    {
        int id1 = this->GetConnPointValue( i );
        int id2 = this->GetConnDonorPointValue( i );

        CgIntField fNodeId1, fNodeId2;
        sZone->GetElementNodeId( id1, fNodeId1 );
        tZone->GetElementNodeId( id2, fNodeId2 );

        f2fmap.AddFacePoint(fNodeId1, fNodeId2, nodeMesh1, nodeMesh2 );
    }
}

#endif

EndNameSpace
