/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "CgnsFile.h"
#include "NodeMesh.h"
#include "Dimension.h"
#include <iostream>
#include <iomanip>

BeginNameSpace( ONEFLOW )
#ifdef ENABLE_CGNS

CgnsCoor::CgnsCoor( CgnsZone * cgnsZone )
{
    this->cgnsZone = cgnsZone;
    this->ndim = 3;
    this->typeList.resize( this->ndim );
    this->coor.resize( this->ndim );
    this->nNodeList.resize( this->ndim );
    this->coorNameList.resize( this->ndim );
    this->nodeMesh = new NodeMesh();
}

CgnsCoor::~CgnsCoor()
{
    DeAlloc();
    delete this->nodeMesh;
}

CgInt CgnsCoor::GetNNode()
{
    return this->nNodes;
}

CgInt CgnsCoor::GetNCell()
{
    return this->nCells;
}

void CgnsCoor::SetNNode( CgInt nNodes )
{
    this->nNodes = nNodes;
}

void CgnsCoor::SetNCell( CgInt nCells )
{
    this->nCells = nCells;
}

void CgnsCoor::Alloc( int iCoor, int nNodes, DataType_t data_type )
{
    if ( data_type == RealSingle )
    {
        this->coor[ iCoor ] = new float [ nNodes ];
    }
    else
    {
        this->coor[ iCoor ] = new double [ nNodes ];
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
    int nNodes = this->nNodeList[ iCoor ];
    if ( data_type == RealSingle )
    {
        float * data = static_cast< float * >( this->coor[ iCoor ] );
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            var[ iNode ] = data[ iNode ];
        }
    }
    else
    {
        double * data = static_cast< double * >( this->coor[ iCoor ] );
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            var[ iNode ] = data[ iNode ];
        }
    }
}

void CgnsCoor::SetCoorData( int iCoor, DataType_t data_type, Real * var )
{
    int nNodes = this->nNodeList[ iCoor ];
    if ( data_type == RealSingle )
    {
        float * data = static_cast< float * >( this->coor[ iCoor ] );
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            data[ iNode ] = var[ iNode ];
        }
    }
    else
    {
        double * data = static_cast< double * >( this->coor[ iCoor ] );
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            data[ iNode ] = var[ iNode ];
            var[ iNode ] = data[ iNode ];
        }
    }
}

void CgnsCoor::SetAllCoorData()
{
    HXVector< Real * > xyz( this->ndim );

    xyz[ 0 ] = & nodeMesh->xN[ 0 ];
    xyz[ 1 ] = & nodeMesh->yN[ 0 ];
    xyz[ 2 ] = & nodeMesh->zN[ 0 ];

    for ( int iCoor = 0; iCoor < this->ndim; ++ iCoor )
    {
        DataType_t data_type = this->typeList[ iCoor ];
        this->SetCoorData( iCoor, data_type, xyz[ iCoor ] );
    }
}

void CgnsCoor::CopyCoorData( CgnsCoor * cgnsCoorIn )
{
    for ( int iCoor = 0; iCoor < this->nCoor; ++ iCoor )
    {
        int nNodes = this->nNodeList[ iCoor ];
        DataType_t dataType = this->typeList[ iCoor ];
        this->Alloc( iCoor, nNodes, dataType );

        if ( dataType == RealSingle )
        {
            float * data = static_cast<float *>( this->coor[ iCoor ] );
            float * dataIn = static_cast<float *>( cgnsCoorIn->coor[ iCoor ] );
            for ( int iNode = 0; iNode < nNodes; ++ iNode )
            {
                data[ iNode ] = dataIn[ iNode ];
            }
        }
        else
        {
            double * data = static_cast<double *>( this->coor[ iCoor ] );
            double * dataIn = static_cast<double *>( cgnsCoorIn->coor[ iCoor ] );
            for ( int iNode = 0; iNode < nNodes; ++ iNode )
            {
                data[ iNode ] = dataIn[ iNode ];
            }
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
    int fileId = this->cgnsZone->cgnsBase->cgnsFile->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zoneId = this->cgnsZone->zId;

    cg_ncoords( fileId, baseId, zoneId, & this->nCoor );
    std::cout << "   this->nCoor = " << this->nCoor << "\n";

    int nNodes = this->GetNNode();

    for ( int iCoor = 0; iCoor < this->nCoor; ++ iCoor )
    {
        DataType_t dataType;
        CgnsTraits::char33 coorName;
        int coordId = iCoor + 1;
        cg_coord_info( fileId, baseId, zoneId, coordId, & dataType, coorName );
        std::cout << "   coorName = " << coorName << " dataType = " << dataType << " dataTypeName = " << DataTypeName[ dataType ] << "\n";
        this->typeList[ iCoor ] = dataType;
        this->nNodeList[ iCoor ] = nNodes;
        this->coorNameList[ iCoor ] = coorName;
        this->Alloc( iCoor, static_cast<int>( nNodes ), dataType );
        //Read the x-, y-, z-coordinates.
        cg_coord_read( fileId, baseId, zoneId, coorName, dataType, this->irmin, this->irmax, this->GetCoor( iCoor ) );
    }

    NodeMesh * nodeMesh = this->GetNodeMesh();
    nodeMesh->CreateNodes( static_cast<int>( nNodes ) );

    this->SetAllData( nodeMesh->xN, nodeMesh->yN, nodeMesh->zN );
}

void CgnsCoor::ReadCgnsGridCoordinates( CgnsCoor * cgnsCoorIn )
{
    //Determine the number and names of the coordinates.
    int fileId = this->cgnsZone->cgnsBase->cgnsFile->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zoneId = this->cgnsZone->zId;

    std::cout << " this->nCoor = " << this->nCoor << "\n";
    this->nCoor = cgnsCoorIn->nCoor;
    std::cout << " this->nCoor = " << this->nCoor << "\n";

    int nNodes = this->GetNNode();

    std::cout << " this->nNodes = " << this->nNodes << "\n";

    this->typeList = cgnsCoorIn->typeList;
    this->nNodeList = cgnsCoorIn->nNodeList;
    this->coorNameList = cgnsCoorIn->coorNameList;
    this->CopyCoorData( cgnsCoorIn );

    NodeMesh * nodeMesh = this->GetNodeMesh();
    NodeMesh * nodeMeshIn = cgnsCoorIn->GetNodeMesh();

    * nodeMesh = * nodeMeshIn;
}

void CgnsCoor::DumpCgnsGridCoordinates()
{
    //Determine the number and names of the coordinates.
    int fileId = this->cgnsZone->cgnsBase->cgnsFile->fileId;
    int baseId = this->cgnsZone->cgnsBase->baseId;
    int zoneId = this->cgnsZone->zId;

     std::cout << "   this->nCoor = " << this->nCoor << "\n";

     for ( int iCoor = 0; iCoor < this->nCoor; ++ iCoor )
     {
         DataType_t dataType = this->typeList[ iCoor ];
         std::string & coorName = this->coorNameList[ iCoor ];
          //Write the x-, y-, z-coordinates.
         int indexCoor = -1;
         cg_coord_write( fileId, baseId, zoneId, dataType, coorName.c_str(), this->GetCoor( iCoor ), &indexCoor );
         std::cout << "   coorName = " << coorName << " dataType = " << dataType << " indexCoor = " << indexCoor << "\n";
     }
}

void CgnsCoor::FreeMesh()
{
    delete this->nodeMesh;
    this->nodeMesh = 0;
}

NodeMesh * CgnsCoor::GetNodeMesh()
{
    return this->nodeMesh;
}

void CgnsCoor::SetDimension()
{
    int rind[ 6 ];
    int result = cg_rind_read( rind );
    if ( result != CG_OK )
    {
        for ( int i = 0; i < 6; ++ i )
        {
            rind[i] = 0;
        }
    }

    CgInt * isize = this->cgnsZone->isize;

    if ( this->cgnsZone->cgnsZoneType == CGNS_ENUMV( Structured ) )
    {
        this->SetDimensionStr();
    }
    else
    {
        this->SetDimensionUns();
    }

    std::cout << "   numberOfNodes = " << this->GetNNode() << " numberOfCells = " << this->GetNCell() << "\n";
}

void CgnsCoor::SetDimensionStr()
{
    CgInt * isize = this->cgnsZone->isize;

    // lower range index
    irmin[ 0 ] = 1;
    irmin[ 1 ] = 1;
    irmin[ 2 ] = 1;

    // upper range index of vertices
    irmax[ 0 ] = 1;
    irmax[ 1 ] = 1;
    irmax[ 2 ] = 1;

    cellSize[ 0 ] = 1;
    cellSize[ 1 ] = 1;
    cellSize[ 2 ] = 1;

    // upper range index of vertices
    // vertex size
    int j = 0;
    irmax[ 0 ] = isize[ j ++ ];
    irmax[ 1 ] = isize[ j ++ ];
    if ( this->cgnsZone->cgnsBase->celldim == THREE_D )
    {
        irmax[ 2 ] = isize[ j ++ ];
    }
    // cell size
    cellSize[ 0 ] = isize[ j ++ ];
    cellSize[ 1 ] = isize[ j ++ ];
    if ( this->cgnsZone->cgnsBase->celldim == THREE_D )
    {
        cellSize[ 2 ] = isize[ j ++ ];
    }
    std::cout << "   The Dimension Of Grid is : \n";
    std::cout << "   I Direction " << std::setw( 10 ) << irmin[ 0 ] << std::setw( 10 ) << irmax[ 0 ] << "\n";
    std::cout << "   J Direction " << std::setw( 10 ) << irmin[ 1 ] << std::setw( 10 ) << irmax[ 1 ] << "\n";
    std::cout << "   K Direction " << std::setw( 10 ) << irmin[ 2 ] << std::setw( 10 ) << irmax[ 2 ] << "\n";
    int nNodes = irmax[ 0 ] * irmax[ 1 ] * irmax[ 2 ];
    int nCells = cellSize[ 0 ] * cellSize[ 1 ] * cellSize[ 2 ];
    this->SetNNode( nNodes );
    this->SetNCell( nCells );
}

void CgnsCoor::SetDimensionUns()
{
    CgInt * isize = this->cgnsZone->isize;
    irmin[ 0 ] = 1;
    irmin[ 1 ] = 0;
    irmin[ 2 ] = 0;

    irmax[ 0 ] = isize[ 0 ];
    irmax[ 1 ] = 0;
    irmax[ 2 ] = 0;

    cellSize[ 0 ] = isize[ 1 ];

    int nNodes = irmax[ 0 ];
    int nCells = cellSize[ 0 ];
    this->SetNNode( nNodes );
    this->SetNCell( nCells );
}

void CgnsCoor::SetDimension( CgnsCoor * cgnsCoorIn )
{
    CgInt * isize = this->cgnsZone->isize;

    isize[ 0 ] = cgnsCoorIn->GetNNode();
    isize[ 1 ] = cgnsCoorIn->GetNCell();

    irmin[ 0 ] = 1;
    irmin[ 1 ] = 0;
    irmin[ 2 ] = 0;

    irmax[ 0 ] = isize[ 0 ];
    irmax[ 1 ] = 0;
    irmax[ 2 ] = 0;

    cellSize[ 0 ] = isize[ 1 ];

    this->SetNNode( irmax[ 0 ] );
    this->SetNCell( cellSize[ 0 ] );

    this->cgnsZone->InitLgMapping();

    std::cout << "   numberOfNodes = " << this->GetNNode() << " numberOfCells = " << this->GetNCell() << "\n";
}


#endif
EndNameSpace
