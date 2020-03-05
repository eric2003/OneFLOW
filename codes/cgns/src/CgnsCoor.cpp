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

#include "CgnsCoor.h"
#include "CgnsZone.h"
#include "CgnsBase.h"
#include "NodeMesh.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsCoor::CgnsCoor( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->ndim = 3;
    this->typeList.resize( this->ndim );
    this->coor.resize( this->ndim );
    this->nNodeList.resize( this->ndim );
    this->nodeMesh = new NodeMesh();
}

CgnsCoor::~CgnsCoor()
{
    DeAlloc();
    delete this->nodeMesh;
}

void CgnsCoor::Alloc( int iCoor, int nNode, DataType_t data_type )
{
    this->typeList[ iCoor ] = data_type;
    this->nNodeList[ iCoor ] = nNode;

    if ( data_type == RealSingle )
    {
        this->coor[ iCoor ] = new float [ nNode ];
    }
    else
    {
        this->coor[ iCoor ] = new double [ nNode ];
    }
}

void CgnsCoor::SetAllData( RealField & x, RealField & y, RealField & z )
{
    HXVector< Real * > xyz( this->ndim );

    xyz[ 0 ] = & x[ 0 ];
    xyz[ 1 ] = & y[ 0 ];
    xyz[ 2 ] = & z[ 0 ];

    for ( int iCoor = 0; iCoor < this->ndim; ++ iCoor )
    {
        DataType_t data_type = this->typeList[ iCoor ];
        SetData( iCoor, data_type, xyz[ iCoor ] );
    }
}

void CgnsCoor::SetData( int iCoor, DataType_t data_type, Real * var )
{
    int nNode = this->nNodeList[ iCoor ];
    if ( data_type == RealSingle )
    {
        float * data = static_cast< float * >( this->coor[ iCoor ] );
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            var[ iNode ] = data[ iNode ];
        }
    }
    else
    {
        double * data = static_cast< double * >( this->coor[ iCoor ] );
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            var[ iNode ] = data[ iNode ];
        }
    }
}

void CgnsCoor::DeAlloc()
{
    for ( int iCoor = 0; iCoor < this->ndim; ++ iCoor )
    {
        int data_type = this->typeList[ iCoor ];
        if ( data_type == RealSingle )
        {
            float * data  = static_cast< float * >( this->coor[ iCoor ] );
            delete [] data;
        }
        else
        {
            double * data = static_cast< double * >( this->coor[ iCoor ] );
            delete [] data;
        }
    }
}

void CgnsCoor::ReadCgnsGridCoordinates()
{
    //Determine the number and names of the coordinates.
    int fileId = this->cgnsZone->cgnsBase->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zoneId = this->cgnsZone->zId;

    cg_ncoords( fileId, baseId, zoneId, & this->nCoor );

    int nNode = this->cgnsZone->GetNNode();

    for ( int coordId = 0; coordId < this->nCoor; ++ coordId )
    {
        DataType_t dataType;
        CgnsTraits::char33 coorName;
        cg_coord_info( fileId, baseId, zoneId, coordId + 1, & dataType, coorName );
        this->Alloc( coordId, static_cast<int>( nNode ), dataType );
        //Read the x-, y-, z-coordinates.
        cg_coord_read( fileId, baseId, zoneId, coorName, dataType, this->cgnsZone->irmin, this->cgnsZone->irmax, this->GetCoor( coordId ) );
    }
}

void CgnsCoor::FreeMesh()
{
    delete this->nodeMesh;
    this->nodeMesh = 0;
}


#endif
EndNameSpace