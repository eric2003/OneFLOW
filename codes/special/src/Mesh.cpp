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

#include "Mesh.h"
#include "HXSort.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "NodeMesh.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "BcRecord.h"
#include "ElementHome.h"
#include "HXMath.h"
#include "HXCgns.h"
#include "Visual.h"
#include "Dimension.h"
#include "DataBase.h"
#include <algorithm>
#include <iostream>
#include <ctime>
using namespace std;

BeginNameSpace( ONEFLOW )

HXRandomClass::HXRandomClass()
{
}

HXRandomClass::~HXRandomClass()
{
}

void HXRandomClass::Initialize()
{
    ::srand( ::time( 0 ) );
}

int HXRandomClass::Random( int rangeMin, int rangeMax )
{
    double value = static_cast< double > ( rand() ) / RAND_MAX * ( rangeMax - rangeMin + 1 ) + rangeMin;
    return static_cast< int > ( value );
}

void HXRandomClass::RangeRandom( int rangeMin, int rangeMax, vector< int > & results )
{
    int numberOfElements = results.size();
    vector< int > flag( numberOfElements, 0 );
    for ( int iElement = 0; iElement< numberOfElements; ++ iElement )
    {
        int candidateNumber = - 1;
        while ( true )
        {
            candidateNumber = HXRandomClass::Random( rangeMin, rangeMax );
            if ( flag[ candidateNumber ] == 0 )
            {
                flag[ candidateNumber ] = 1;
                results[ iElement ] = candidateNumber;
                break;
            }
        }
    }
}


SimpleMesh2D::SimpleMesh2D()
{
    ;
}

SimpleMesh2D::~SimpleMesh2D()
{
    ;
}

void SimpleMesh2D::GenerateMesh()
{
    this->GenerateCircleMesh();
}

void SimpleMesh2D::GenerateCircleMesh()
{
    int ni = 36;
    int nj = 10;

    RealField x1Array( ni ), y1Array( ni );
    RealField x2Array( ni ), y2Array( ni );
    IntField nodeArray1( ni );
    IntField nodeArray2( ni );

    this->GenerateCircleSurface( x1Array, y1Array, ni );

    this->PushCircleNode( x1Array, y1Array, nodeArray1 );

    for ( int j = 0; j < nj; ++ j )
    {
        this->CalcX2Y2Array( x1Array, y1Array, x2Array, y2Array );

        this->PushCircleNode( x2Array, y2Array, nodeArray2 );

        this->PushElement( nodeArray1, nodeArray2 );

        x1Array = x2Array;
        y1Array = y2Array;

        nodeArray1 = nodeArray2;
    }
}

void SimpleMesh2D::GenerateExtrapolationMesh()
{
    int ni = 11;
    int nj = 1;

    RealField x1Array( ni ), y1Array( ni );
    RealField x2Array( ni ), y2Array( ni );
    IntField nodeArray1( ni );
    IntField nodeArray2( ni );

    for ( int i = 0; i < ni; ++ i )
    {
        x1Array[ i ] = i;
        y1Array[ i ] = 0;
    }

    this->PushCircleNode( x1Array, y1Array, nodeArray1 );

    for ( int j = 0; j < nj; ++ j )
    {
        this->CalcX2Y2ArrayNoLoop( x1Array, y1Array, x2Array, y2Array );
        this->PushCircleNode( x2Array, y2Array, nodeArray2 );
        this->PushElement( nodeArray1, nodeArray2, 1 );

        x1Array = x2Array;
        y1Array = y2Array;
        nodeArray1 = nodeArray2;
    }
}

void SimpleMesh2D::GenerateRectangleMesh()
{
    //construct rectangle region
    this->ni = 21;
    this->nj = 21;

    ONEFLOW::AllocateVector( xx, ni, nj );
    ONEFLOW::AllocateVector( yy, ni, nj );
    ONEFLOW::AllocateVector( zz, ni, nj );

    //Real xmin = - 4.0;
    //Real xmax =   4.0;
    //Real ymin = - 4.0;
    //Real ymax =   4.0;

    Real xmin = - 10.0;
    Real xmax = 10.0;
    Real ymin = - 10.0;
    Real ymax = 10.0;


    Real dx = ( xmax - xmin ) / ( ni - 1 );
    Real dy = ( ymax - ymin ) / ( nj - 1 );

    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            xx[ i ][ j ] = xmin + i * dx;
            yy[ i ][ j ] = ymin + j * dy;
            zz[ i ][ j ] = 0.0;
        }
    }

    ONEFLOW::AllocateVector( ijkNodeMapping, ni, nj );

    int iCount = 0;
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            ijkNodeMapping[ i ][ j ] = iCount ++;
            Real xNode = this->xx[ i ][ j ];
            Real yNode = this->yy[ i ][ j ];
            Real zNode = this->zz[ i ][ j ];
            mesh->nodeMesh->xN.push_back( xNode );
            mesh->nodeMesh->yN.push_back( yNode );
            mesh->nodeMesh->zN.push_back( zNode );
        }
    }

    this->ConstructElement();
}

void SimpleMesh2D::ConstructElement()
{
    HXRandomClass::Initialize();
    HXRandomClass::Random( 0, 1 );

    for ( int j = 0; j < nj - 1; ++ j )
    {
        for ( int i = 0; i < ni - 1; ++ i )
        {
            int p1 = this->ijkNodeMapping[ i ][ j ];
            int p2 = this->ijkNodeMapping[ i + 1 ][ j ];
            int p3 = this->ijkNodeMapping[ i + 1 ][ j + 1 ];
            int p4 = this->ijkNodeMapping[ i ][ j + 1 ];

            int typeFlag = HXRandomClass::Random( 0, 1 );

            int elementType = ONEFLOW::QUAD_4;

            mesh->cellMesh->cellTopo->PushElement( p1, p2, p3, p4, elementType );
        }
    }
}

void SimpleMesh2D::PushElement( IntField & nodeArray1, IntField & nodeArray2, int shift )
{
    UInt numberOfNodes = nodeArray1.size();
    for ( UInt iNode = 0; iNode < numberOfNodes - shift; ++ iNode )
    {
        int iNode1 = iNode;
        int iNode2 = ( iNode + 1 ) % numberOfNodes;
        int p1 = nodeArray1[ iNode1 ];
        int p2 = nodeArray1[ iNode2 ];
        int p3 = nodeArray2[ iNode1 ];
        int p4 = nodeArray2[ iNode2 ];
        this->mesh->cellMesh->cellTopo->PushElement( p1, p2, p3, ONEFLOW::TRI_3 );
        this->mesh->cellMesh->cellTopo->PushElement( p2, p4, p3, ONEFLOW::TRI_3 );
    }
}

void SimpleMesh2D::GenerateCircleSurface( RealField & xArray, RealField & yArray, int ni )
{
    Real dcit = 2.0 * ONEFLOW::PI / ni;
    Real r0 = 1.0;
    Real cit0 = 0.0;
    //The premise is clockwise
    for ( int i = 0; i < ni; ++ i )
    {
        Real cit = cit0 - i * dcit;
        Real x = r0 * cos( cit );
        Real y = r0 * sin( cit );
        xArray[ i ] = x;
        yArray[ i ] = y;
    }
}

void SimpleMesh2D::CalcX2Y2Array( RealField & x1Array, RealField & y1Array, RealField & x2Array, RealField & y2Array )
{
    int numberOfNodes = x1Array.size();
    Real factor = sin( PI / 3 );
    for ( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        int iNode1 = iNode;
        int iNode2 = ( iNode + 1 ) % numberOfNodes;
        Real x1 = x1Array[ iNode1 ];
        Real y1 = y1Array[ iNode1 ];
        Real x2 = x1Array[ iNode2 ];
        Real y2 = y1Array[ iNode2 ];
        Real xm = 0.5 * ( x1 + x2 );
        Real ym = 0.5 * ( y1 + y2 );
        Real dx = x2 - x1;
        Real dy = y2 - y1;
        Real ds = ONEFLOW::DIST( dx, dy );
        Real nx = dy / ds;
        Real ny = - dx / ds;

        Real dr = ds * factor;

        Real x = xm - nx * dr;
        Real y = ym - ny * dr;

        x2Array[ iNode ] = x;
        y2Array[ iNode ] = y;
    }
}

void SimpleMesh2D::CalcX2Y2ArrayNoLoop( RealField & x1Array, RealField & y1Array, RealField & x2Array, RealField & y2Array )
{
    int numberOfNodes = x1Array.size();
    Real factor = sin( PI / 3 );
    for ( int iNode = 0; iNode < numberOfNodes - 1; ++ iNode )
    {
        int iNode1 = iNode;
        int iNode2 = ( iNode + 1 ) % numberOfNodes;
        Real x1 = x1Array[ iNode1 ];
        Real y1 = y1Array[ iNode1 ];
        Real x2 = x1Array[ iNode2 ];
        Real y2 = y1Array[ iNode2 ];
        Real xm = 0.5 * ( x1 + x2 );
        Real ym = 0.5 * ( y1 + y2 );

        Real dx = x2 - x1;
        Real dy = y2 - y1;
        Real ds = ONEFLOW::DIST( dx, dy );
        Real nx = dy / ds;
        Real ny = - dx / ds;

        Real dr = ds * factor;

        Real x = xm - nx * dr;
        Real y = ym - ny * dr;
        x2Array[ iNode ] = x;
        y2Array[ iNode ] = y;
    }

    int nNode = numberOfNodes - 1;

    x2Array[ nNode ] = 2 * x2Array[ nNode - 1 ] - x2Array[ nNode - 2 ];
    y2Array[ nNode ] = 2 * y2Array[ nNode - 1 ] - y2Array[ nNode - 2 ];
}

void SimpleMesh2D::PushCircleNode( RealField & xArray, RealField & yArray, IntField & nodeArray )
{
    int numberOfNodes = xArray.size();
    for ( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        nodeArray[ iNode ] = mesh->nodeMesh->xN.size();
        mesh->nodeMesh->xN.push_back( xArray[ iNode ] );
        mesh->nodeMesh->yN.push_back( yArray[ iNode ] );
        mesh->nodeMesh->zN.push_back( 0.0 );
    }
}

Mesh::Mesh()
{
    nodeMesh = 0;
    faceMesh = 0;
    cellMesh = 0;
    dataBase = new DataBase();
}

Mesh::~Mesh()
{
    delete nodeMesh;
    delete faceMesh;
    delete cellMesh;
    delete dataBase;
}

void Mesh::CreateMesh()
{
    cout << "Mesh::CreateMesh()\n";
    nodeMesh = new NodeMesh();
    faceMesh = new FaceMesh();
    cellMesh = new CellMesh();

    SimpleMesh2D simpleMesh2D;
    simpleMesh2D.SetMesh( this );
    simpleMesh2D.GenerateMesh();
    this->ConstructTopology();
    ONEFLOW::Visual::Show( this );

    this->CalcMetrics();
}

void Mesh::ConstructTopology()
{
    UInt numberOfNodes = this->nodeMesh->GetNumberOfNodes();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();
    set< HXSort< IntField > > faceSet;
    HXSort< IntField > faceForSorting;

    CellTopo * cellTopo = this->cellMesh->cellTopo;
    FaceTopo * faceTopo = this->faceMesh->faceTopo;

    for ( UInt iCell = 0; iCell < numberOfCells; ++ iCell )
    {
        IntField & element = cellTopo->cellToNode[ iCell ];

        int elementType = cellTopo->cellType[ iCell ];

        UnitElement * unitElement = ONEFLOW::ElementHome::GetUnitElement( elementType );
        int numberOfFaceInElement = unitElement->GetElementFaceNumber();

        for ( int iLocalFace = 0; iLocalFace < numberOfFaceInElement; ++ iLocalFace )
        {
            IntField & localFaceNodeIndexArray = unitElement->GetElementFace( iLocalFace );
            int faceType = unitElement->GetFaceType( iLocalFace );
            int numberOfFacePoints = localFaceNodeIndexArray.size();
            IntField faceNodeIndexArray;
            for ( int iFacePoint = 0; iFacePoint < numberOfFacePoints; ++ iFacePoint )
            {
                faceNodeIndexArray.push_back( element[ localFaceNodeIndexArray[ iFacePoint ] ] );
            }

            IntField faceNodeIndexArraySort = faceNodeIndexArray;
            std::sort( faceNodeIndexArraySort.begin(), faceNodeIndexArraySort.end() );
            faceForSorting.value = faceNodeIndexArraySort;

            set< HXSort< IntField > >::iterator iter = faceSet.find( faceForSorting );
            if ( iter == faceSet.end() )
            {
                faceForSorting.index = faceSet.size();
                faceSet.insert( faceForSorting );
                int faceIndex = faceForSorting.index;
                int newSize = faceIndex + 1;
                faceTopo->lCell.resize( newSize );
                faceTopo->rCell.resize( newSize );
                faceTopo->lPosition.resize( newSize );
                faceTopo->rPosition.resize( newSize );
                faceTopo->faceType.resize( newSize );

                faceTopo->faceType[ faceIndex ] = faceType;
                faceTopo->lCell[ faceIndex ] = iCell;
                faceTopo->rCell[ faceIndex ] = ONEFLOW::INVALID_INDEX;
                faceTopo->lPosition[ faceIndex ] = iLocalFace;
                faceTopo->rPosition[ faceIndex ] = ONEFLOW::INVALID_INDEX;
                faceTopo->f2n.push_back( faceNodeIndexArray );
            }
            else
            {
                int faceIndex = iter->index;
                faceTopo->rCell[ faceIndex ] = iCell;
                faceTopo->rPosition[ faceIndex ] = iLocalFace;
            }
        }
    }

    this->SwapBoundary();
}

void Mesh::SwapBoundary()
{
    UInt nFace = this->faceMesh->GetNFace();

    IntField orderMapping( nFace );

    CellTopo * cellTopo = this->cellMesh->cellTopo;
    FaceTopo * faceTopo = this->faceMesh->faceTopo;

    int iBoundaryFaceCount = 0;
    int iCount = 0;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int rc = faceTopo->rCell[ iFace ];
        if ( rc == ONEFLOW::INVALID_INDEX )
        {
            orderMapping[ iCount ++ ] = iFace;
            ++ iBoundaryFaceCount;
        }
    }

    int nBFace = iBoundaryFaceCount;
    this->faceMesh->SetNBFace( nBFace );

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int rc = faceTopo->rCell[ iFace ];
        if ( rc != ONEFLOW::INVALID_INDEX )
        {
            orderMapping[ iCount ++ ] = iFace;
        }
    }

    IntField lCellIndexSwap = faceTopo->lCell;
    IntField rCellIndexSwap = faceTopo->rCell;

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int oldFaceIndex = orderMapping[ iFace ];
        faceTopo->lCell[ iFace ] = lCellIndexSwap[ oldFaceIndex ];
        faceTopo->rCell[ iFace ] = rCellIndexSwap[ oldFaceIndex ];
    }

    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        faceTopo->rCell[ iFace ] = iFace + numberOfCells;
    }

    IntField lPositionSwap = faceTopo->lPosition;
    IntField rPositionSwap = faceTopo->rPosition;

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int oldFaceIndex = orderMapping[ iFace ];
        faceTopo->lPosition[ iFace ] = lPositionSwap[ oldFaceIndex ];
        faceTopo->rPosition[ iFace ] = rPositionSwap[ oldFaceIndex ];
    }

    LinkField faceToNodeSwap = faceTopo->f2n;

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int oldFaceIndex = orderMapping[ iFace ];
        faceTopo->f2n[ iFace ] = faceToNodeSwap[ oldFaceIndex ];
    }

    IntField faceTypeSwap = faceTopo->faceType;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int oldFaceIndex = orderMapping[ iFace ];
        faceTopo->faceType[ iFace ] = faceTypeSwap[ oldFaceIndex ];
    }
}

void Mesh::AllocateMetrics()
{
    this->faceMesh->AllocateMetrics();
    this->cellMesh->AllocateMetrics( this->faceMesh );
}

void Mesh::CalcMetrics()
{
    this->AllocateMetrics();
    if ( ONEFLOW::IsOneD() )
    {
        this->CalcMetrics1D();
    }
    else if ( ONEFLOW::IsTwoD() )
    {
        this->CalcMetrics2D();
    }
    else if ( ONEFLOW::IsThreeD() )
    {
        this->CalcMetrics3D();
    }
}

void Mesh::CalcMetrics1D()
{
    //compute face center first for one dimensional case then face normal
    this->CalcFaceCenter1D();
    this->CalcCellCenterVol1D();
    this->CalcFaceNormal1D();
    this->CalcGhostCellCenterVol1D();
}

void Mesh::CalcMetrics2D()
{
    this->CalcFaceNormal2D();
    this->CalcFaceCenter2D();
    this->CalcCellCenterVol2D();
}

void Mesh::CalcMetrics3D()
{
    this->CalcFaceNormal3D();
    this->CalcFaceCenter3D();
    this->CalcCellCenterVol3D();
}

void Mesh::CalcFaceNormal2D()
{
    this->faceMesh->CalcFaceNormal2D( this->nodeMesh );
}

void Mesh::CalcFaceCenter2D()
{
    this->faceMesh->CalcFaceCenter2D( this->nodeMesh );
}

void Mesh::CalcFaceCenter1D()
{
    this->faceMesh->CalcFaceCenter1D( this->nodeMesh );
}

void Mesh::CalcFaceNormal1D()
{
    this->faceMesh->CalcFaceNormal1D( this->nodeMesh, this->cellMesh );
}

void Mesh::CalcCellCenterVol1D()
{
    UInt nFace = this->faceMesh->GetNFace();
    UInt nBFace = this->faceMesh->GetNBFace();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    RealField & xcc  = this->cellMesh->xcc ;
    RealField & ycc  = this->cellMesh->ycc ;
    RealField & zcc  = this->cellMesh->zcc ;
    RealField & vol = this->cellMesh->vol;

    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    CellTopo * cellTopo = this->cellMesh->cellTopo;
    FaceTopo * faceTopo = this->faceMesh->faceTopo;

    for ( UInt iCell = 0; iCell < numberOfCells; ++ iCell )
    {
        IntField & element = cellTopo->cellToNode[ iCell ];
        int p1 = element[ 0 ];
        int p2 = element[ 1 ];
        xcc[ iCell  ] = half * ( xN[ p1 ] + xN[ p2 ] );
        ycc[ iCell  ] = half * ( yN[ p1 ] + yN[ p2 ] );
        zcc[ iCell  ] = half * ( zN[ p1 ] + zN[ p2 ] );
        Real dx = xN[ p2 ] - xN[ p1 ];
        Real dy = yN[ p2 ] - yN[ p1 ];
        Real dz = zN[ p2 ] - zN[ p1 ];
        vol[ iCell  ] = ONEFLOW::DIST( dx, dy, dz );
    }
}

void Mesh::CalcGhostCellCenterVol1D()
{
    UInt nFace = this->faceMesh->GetNFace();
    UInt nBFace = this->faceMesh->GetNBFace();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    RealField & xcc  = this->cellMesh->xcc ;
    RealField & ycc  = this->cellMesh->ycc ;
    RealField & zcc  = this->cellMesh->zcc ;
    RealField & vol = this->cellMesh->vol;

    RealField & xfn = this->faceMesh->xfn;
    RealField & yfn = this->faceMesh->yfn;
    RealField & zfn = this->faceMesh->zfn;

    RealField & xfc = this->faceMesh->xfc;
    RealField & yfc = this->faceMesh->yfc;
    RealField & zfc = this->faceMesh->zfc;

    RealField & area = this->faceMesh->area;

    CellTopo * cellTopo = this->cellMesh->cellTopo;
    FaceTopo * faceTopo = this->faceMesh->faceTopo;

    // For ghost cells
    for ( UInt iFace = 0; iFace < nBFace; ++ iFace )
    {
        int lc  = faceTopo->lCell[ iFace ];
        int rc = iFace + numberOfCells;
        if ( area[ iFace ] > SMALL )
        {
            Real tmp = 2.0 * ( ( xcc[ lc ] - xfc[ iFace ] ) * xfn[ iFace ]
                             + ( ycc[ lc ] - yfc[ iFace ] ) * yfn[ iFace ]
                             + ( zcc[ lc ] - zfc[ iFace ] ) * zfn[ iFace ] );
            xcc[ rc ] = xcc[ lc ] - xfn[ iFace ] * tmp;
            ycc[ rc ] = ycc[ lc ] - yfn[ iFace ] * tmp;
            zcc[ rc ] = zcc[ lc ] - zfn[ iFace ] * tmp;
        }
        else
        {
            // Degenerated faces
            xcc[ rc ] = - xcc[ lc ] + 2.0 * xfc[ iFace ];
            ycc[ rc ] = - ycc[ lc ] + 2.0 * yfc[ iFace ];
            zcc[ rc ] = - zcc[ lc ] + 2.0 * zfc[ iFace ];
        }
        vol[ rc ] = vol[ lc ];
    }
}

void Mesh::CalcCellCenterVol2D()
{
    UInt nFace = this->faceMesh->GetNFace();
    UInt nBFace = this->faceMesh->GetNBFace();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    RealField & xcc  = this->cellMesh->xcc ;
    RealField & ycc  = this->cellMesh->ycc ;
    RealField & zcc  = this->cellMesh->zcc ;
    RealField & vol = this->cellMesh->vol;

    RealField & xfn = this->faceMesh->xfn;
    RealField & yfn = this->faceMesh->yfn;
    RealField & zfn = this->faceMesh->zfn;

    RealField & xfc = this->faceMesh->xfc;
    RealField & yfc = this->faceMesh->yfc;
    RealField & zfc = this->faceMesh->zfc;

    RealField & area = this->faceMesh->area;

    CellTopo * cellTopo = this->cellMesh->cellTopo;
    FaceTopo * faceTopo = this->faceMesh->faceTopo;

    xcc  = 0;
    ycc  = 0;
    zcc  = 0;
    vol = 0;

    for ( UInt iFace = 0; iFace < nBFace; ++ iFace )
    {
        int lc = faceTopo->lCell[ iFace ];
        Real dot = ( xfc[ iFace ] * xfn[ iFace ] +
                     yfc[ iFace ] * yfn[ iFace ] +
                     zfc[ iFace ] * zfn[ iFace ] ) * area[ iFace ];
        xcc [ lc ] += xfc[ iFace ] * dot;
        ycc [ lc ] += yfc[ iFace ] * dot;
        zcc [ lc ] += zfc[ iFace ] * dot;
        vol[ lc ] += dot;
    }

    // For interior cell faces
    for ( UInt iFace = nBFace; iFace < nFace; ++ iFace )
    {
        int lc = faceTopo->lCell[ iFace ];
        int rc = faceTopo->rCell[ iFace ];
        Real dot = ( xfc[ iFace ] * xfn[ iFace ] +
                     yfc[ iFace ] * yfn[ iFace ] +
                     zfc[ iFace ] * zfn[ iFace ] ) * area[ iFace ];
        xcc [ lc ] += xfc[ iFace ] * dot;
        ycc [ lc ] += yfc[ iFace ] * dot;
        zcc [ lc ] += zfc[ iFace ] * dot;
        vol[ lc ] += dot;

        xcc [ rc ] -= xfc[ iFace ] * dot;
        ycc [ rc ] -= yfc[ iFace ] * dot;
        zcc [ rc ] -= zfc[ iFace ] * dot;
        vol[ rc ] -= dot;
    }

    UInt numberOfCellsHaveNegativeVolumes = 0;
    Real minvol = LARGE, maxvol = 0.0;
    UInt indexMinv = 0, indexMaxv = 0;
    for ( UInt iCell = 0; iCell < numberOfCells; ++ iCell )
    {
        Real tmp = 1.0 / ( 1.5 * vol[ iCell ] + SMALL );
        xcc [ iCell ] *= tmp;
        ycc [ iCell ] *= tmp;
        zcc [ iCell ] *= tmp;
        vol[ iCell ] *= half;

        if ( minvol > vol[ iCell ] )
        {
            minvol = vol[ iCell ];
            indexMinv = iCell;
        }
        if ( maxvol < vol[ iCell ] )
        {
            maxvol = vol[ iCell ];
            indexMaxv = iCell;
        }
        if ( vol[ iCell ] <= 0.0 )
        {
            vol[ iCell ] = - vol[ iCell ];
            ++ numberOfCellsHaveNegativeVolumes;
        }
    }

    // For ghost cells
    for ( UInt iFace = 0; iFace < nBFace; ++ iFace )
    {
        int lc = faceTopo->lCell[ iFace ];
        int rc = iFace + numberOfCells;
        if ( area[ iFace ] > SMALL )
        {
            Real tmp = 2.0 * ( ( xcc[ lc ] - xfc[ iFace ] ) * xfn[ iFace ]
                             + ( ycc[ lc ] - yfc[ iFace ] ) * yfn[ iFace ]
                             + ( zcc[ lc ] - zfc[ iFace ] ) * zfn[ iFace ] );  
            xcc[ iFace + numberOfCells ] = xcc[ lc ] - xfn[ iFace ] * tmp;
            ycc[ iFace + numberOfCells ] = ycc[ lc ] - yfn[ iFace ] * tmp;
            zcc[ iFace + numberOfCells ] = zcc[ lc ] - zfn[ iFace ] * tmp;
        }
        else
        {
            // Degenerated faces
            xcc[ iFace + numberOfCells ] = - xcc[ lc ] + 2.0 * xfc[ iFace ];
            ycc[ iFace + numberOfCells ] = - ycc[ lc ] + 2.0 * yfc[ iFace ];
            zcc[ iFace + numberOfCells ] = - zcc[ lc ] + 2.0 * zfc[ iFace ];
        }
        vol[ rc ] = vol[ lc ];
    }
}

void Mesh::CalcCellCenterVol3D()
{
    UInt nFace = this->faceMesh->GetNFace();
    UInt nBFace = this->faceMesh->GetNBFace();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    RealField & xcc  = this->cellMesh->xcc ;
    RealField & ycc  = this->cellMesh->ycc ;
    RealField & zcc  = this->cellMesh->zcc ;
    RealField & vol = this->cellMesh->vol;

    RealField & xfn = this->faceMesh->xfn;
    RealField & yfn = this->faceMesh->yfn;
    RealField & zfn = this->faceMesh->zfn;

    RealField & xfc = this->faceMesh->xfc;
    RealField & yfc = this->faceMesh->yfc;
    RealField & zfc = this->faceMesh->zfc;

    RealField & area = this->faceMesh->area;

    CellTopo * cellTopo = this->cellMesh->cellTopo;
    FaceTopo * faceTopo = this->faceMesh->faceTopo;

    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    xcc  = 0;
    ycc  = 0;
    zcc  = 0;
    vol = 0;

    for ( UInt iFace = 0; iFace < nFace; ++ iFace )
    {
        int lc = faceTopo->lCell[ iFace ];
        int rc = faceTopo->rCell[ iFace ];

        IntField & faceIndex = faceTopo->f2n[ iFace ];

        UInt faceNodeNumber = faceIndex.size();
        for ( UInt iNode = 0; iNode < faceNodeNumber; ++ iNode )
        {
            int index1 = iNode;
            int index2 = ( iNode + 1 ) % faceNodeNumber;
            int p2 = faceIndex[ index1 ];
            int p3 = faceIndex[ index2 ];

            Real x21 = xN[ p2 ] - xfc[ iFace ];
            Real y21 = yN[ p2 ] - yfc[ iFace ];
            Real z21 = zN[ p2 ] - zfc[ iFace ];

            Real x31 = xN[ p3 ] - xfc[ iFace ];
            Real y31 = yN[ p3 ] - yfc[ iFace ];
            Real z31 = zN[ p3 ] - zfc[ iFace ];

            Real nx  = y21 * z31 - y31 * z21;
            Real ny  = z21 * x31 - z31 * x21;
            Real nz  = x21 * y31 - x31 * y21;

            Real trix  = xfc[ iFace ] + xN[ p2 ] + xN[ p3 ];
            Real triy  = yfc[ iFace ] + yN[ p2 ] + yN[ p3 ];
            Real triz  = zfc[ iFace ] + zN[ p2 ] + zN[ p3 ];

            // Cell Center and Volume
            Real tmp = nx * trix + ny * triy + nz * triz;
            trix *= tmp;
            triy *= tmp;
            triz *= tmp;

            xcc [ lc  ] += trix;
            ycc [ lc  ] += triy;
            zcc [ lc  ] += triz;
            vol[ lc  ] += tmp;

            xcc [ rc ] -= trix;
            ycc [ rc ] -= triy;
            zcc [ rc ] -= triz;
            vol[ rc ] -= tmp;
        }
    }

    UInt cell = 0;
    Real minvol = LARGE, maxvol = 0.0;
    int indexMinv = 0, indexMaxv = 0;
    for ( UInt iCell = 0; iCell < numberOfCells; ++ iCell )
    {
        Real tmp     = 1.0 / ( 4.0 * vol[ iCell ] + SMALL );
        xcc [ iCell ] *= tmp;
        ycc [ iCell ] *= tmp;
        zcc [ iCell ] *= tmp;
        vol[ iCell ] /= 18.0;

        if ( minvol > vol[ iCell ] )
        {
            minvol     = vol[ iCell ];
            indexMinv = iCell;
        }
        if ( maxvol < vol[ iCell ] )
        {
            maxvol     = vol[ iCell ];
            indexMaxv = iCell;
        }

        if ( vol[ iCell ] <= 0.0 )
        {
            vol[ iCell ] = - vol[ iCell ];
            ++ cell;
        }
    }

    if ( cell ) cout << cell << " cells have negative vols \n";

    // For ghost cells
    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int lc = faceTopo->lCell[ iFace ];
        int rc = iFace + numberOfCells;

        if ( area[ iFace ] > SMALL )
        {
           Real tmp = 2.0 * ( ( xcc[ lc  ] - xfc[ iFace ] ) * xfn[ iFace ] 
                            + ( ycc[ lc  ] - yfc[ iFace ] ) * yfn[ iFace ]
                            + ( zcc[ lc  ] - zfc[ iFace ] ) * zfn[ iFace ] );
           xcc[ rc ] = xcc[ lc  ] - xfn[ iFace ] * tmp;
           ycc[ rc ] = ycc[ lc  ] - yfn[ iFace ] * tmp;
           zcc[ rc ] = zcc[ lc  ] - zfn[ iFace ] * tmp;

        }
        else
        {
            // Degenerated faces
           xcc[ rc ] = - xcc[ lc  ] + 2.0 * xfc[ iFace ];
           ycc[ rc ] = - ycc[ lc  ] + 2.0 * yfc[ iFace ];
           zcc[ rc ] = - zcc[ lc  ] + 2.0 * zfc[ iFace ];
       }
       vol[ rc ] = vol[ lc  ];
    }
}

void Mesh::CalcFaceNormal3D()
{
    this->faceMesh->CalcFaceNormal3D( this->nodeMesh );
}

void Mesh::CalcFaceCenter3D()
{
    this->faceMesh->CalcFaceCenter3D( this->nodeMesh );
}


EndNameSpace