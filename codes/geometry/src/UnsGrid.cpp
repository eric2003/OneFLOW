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

#include "UnsGrid.h"
#include "BcRecord.h"
#include "InterFace.h"
#include "HXMath.h"
#include "HXMathExt.h"
#include "NodeMesh.h"
#include "FaceTopo.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "Mesh.h"
#include "DataBaseIO.h"
#include "DataBook.h"
#include "VirtualFile.h"
#include "IFaceLink.h"
#include "Boundary.h"
#include "Tolerence.h"
#include "Dimension.h"
#include "LogFile.h"
#include <iostream>
#include <algorithm>
#include <iterator>
using namespace std;

BeginNameSpace( ONEFLOW )

REGISTER_GRID( UnsGrid )

UnsGrid * UnsGridCast( Grid * gridIn )
{
    return static_cast< UnsGrid * >( gridIn );
}

UnsGrid::UnsGrid()
{
    this->faceTopo = 0;
    this->faceMesh = 0;
    this->cellMesh = 0;
}

UnsGrid::~UnsGrid()
{
    delete this->faceTopo;
    delete this->faceMesh;
    delete this->cellMesh;
}

void UnsGrid::Init()
{
    this->BasicInit();
    this->faceTopo = new FaceTopo();
    this->faceMesh = new FaceMesh();
    this->cellMesh = new CellMesh();
    faceTopo->grid = this;
    this->faceMesh->faceTopo = this->faceTopo;
}

void UnsGrid::Decode( DataBook * databook )
{
    this->ReadGrid( databook );
}

void UnsGrid::Encode( DataBook * databook )
{
    this->WriteGrid( databook );
}

void UnsGrid::ReadGrid( DataBook * databook )
{
    cout << "Reading unstructured grid data files......\n";
    //Read the number of nodes, number of elements surface and number of elements

    ONEFLOW::HXRead( databook, this->nNode );
    ONEFLOW::HXRead( databook, this->nFace );
    ONEFLOW::HXRead( databook, this->nCell );

    cout << "Grid dimension = " << Dim::dimension << endl;

    cout << " number of nodes    : " << this->nNode << endl;
    cout << " number of surfaces : " << this->nFace << endl;
    cout << " number of elements : " << this->nCell << endl;

    this->nodeMesh->CreateNodes( this->nNode );
    this->cellMesh->cellTopo->Alloc( this->nCell );

    ONEFLOW::HXRead( databook, this->nodeMesh->xN );
    ONEFLOW::HXRead( databook, this->nodeMesh->yN );
    ONEFLOW::HXRead( databook, this->nodeMesh->zN );

    cout << "The grid nodes have been read\n";
    ONEFLOW::HXRead( databook, this->volBcType  );

    this->nodeMesh->CalcMinMaxBox();
    this->ReadGridFaceTopology( databook );
    this->ReadBoundaryTopology( databook );
    this->NormalizeBc();

    cout << "All the computing information is ready!\n";
}

void UnsGrid::NormalizeBc()
{
    for ( int iFace = 0; iFace < this->nBFace; ++ iFace )
    {
        this->faceTopo->rCell[ iFace ] = iFace + this->nCell;
    }
}

void UnsGrid::ReadGridFaceTopology( DataBook * databook )
{
    this->faceTopo->f2n.resize( this->nFace );
    this->faceTopo->lCell.resize( this->nFace );
    this->faceTopo->rCell.resize( this->nFace );
    this->faceTopo->faceType.resize( this->nFace );

    IntField numFaceNode( this->nFace );

    ONEFLOW::HXRead( databook, numFaceNode );

    int nsum = ONEFLOW::SUM( numFaceNode );

    cout << "Setting the connection mode of face to point......\n";
    IntField faceNodeMem( nsum );

    ONEFLOW::HXRead( databook, faceNodeMem );

    int ipos = 0;
    for ( int iFace = 0; iFace < this->nFace; ++ iFace )
    {
        int nNode = numFaceNode[ iFace ];
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            int pid = faceNodeMem[ ipos ++ ];
            this->faceTopo->f2n[ iFace ].push_back( pid );
        }
    }

    cout << "Setting the connection mode of face to cell......\n";

    ONEFLOW::HXRead( databook, this->faceTopo->lCell );
    ONEFLOW::HXRead( databook, this->faceTopo->rCell );

    for ( int iFace = 0; iFace < this->nFace; ++ iFace )
    {
        if ( this->faceTopo->lCell[ iFace ] < 0 )
        {
            //need to reverse the node ordering
            IntField & f2n = this->faceTopo->f2n[ iFace ];
            //std::reverse( std::begin( f2n ), std::end( f2n ) );
            std::reverse( f2n.begin(), f2n.end() );
            // now reverse leftCellIndex  and rightCellIndex
            ONEFLOW::SWAP( this->faceTopo->lCell[ iFace ], this->faceTopo->rCell[ iFace ] );
        }
    }
}

void UnsGrid::ReadBoundaryTopology( DataBook * databook )
{
    cout << "Setting the boundary condition......\n";
    ONEFLOW::HXRead( databook, this->nBFace );
    this->faceTopo->SetNBFace( this->nBFace );

    //cout << " nBFace = " << this->nBFace << endl;

    //设置边界条件
    BcRecord * bcRecord = this->faceTopo->bcManager->bcRecord;
    ONEFLOW::HXRead( databook, bcRecord->bcType );
    ONEFLOW::HXRead( databook, bcRecord->bcRegion );
    ONEFLOW::HXRead( databook, this->nIFace );
    cout << " nBFace = " << this->nBFace;
    cout << " nIFace = " << this->nIFace << endl;
    this->interFace->Set( this->nIFace, this );

    if ( this->nIFace > 0 )
    {
        ONEFLOW::HXRead( databook, this->interFace->zoneId            );
        ONEFLOW::HXRead( databook, this->interFace->localInterfaceId  );
        ONEFLOW::HXRead( databook, this->interFace->i2b               );
    }
}

void UnsGrid::WriteBoundaryTopology( VirtualFile * vf )
{
    cout << "Output boundary condition\n";
    BcRecord * bcRecord = this->faceTopo->bcManager->bcRecord;
    int nBFace = bcRecord->GetNBFace();
    cout << "nBFace = " << nBFace << " nIFace = " << this->nIFace << endl;
    cout << "this->interFace->nIFace = " << this->interFace->nIFace << endl;

    ONEFLOW::HXWrite( vf, nBFace );
    ONEFLOW::HXWrite( vf, bcRecord->bcType );
    ONEFLOW::HXWrite( vf, bcRecord->bcRegion );
    ONEFLOW::HXWrite( vf, this->interFace->nIFace );

    if ( this->interFace->nIFace > 0 )
    {
        ONEFLOW::HXWrite( vf, this->interFace->zoneId );
        ONEFLOW::HXWrite( vf, this->interFace->localInterfaceId );
        ONEFLOW::HXWrite( vf, this->interFace->i2b );
    }
}

void UnsGrid::WriteGrid( DataBook * databook )
{
    ONEFLOW::HXWrite( databook, this->nNode );
    ONEFLOW::HXWrite( databook, this->nFace );
    ONEFLOW::HXWrite( databook, this->nCell );

    ONEFLOW::HXWrite( databook, this->nodeMesh->xN );
    ONEFLOW::HXWrite( databook, this->nodeMesh->yN );
    ONEFLOW::HXWrite( databook, this->nodeMesh->zN );

    ONEFLOW::HXWrite( databook, this->volBcType  );

    this->WriteGridFaceTopology( databook );
    this->WriteBoundaryTopology( databook );
}

void UnsGrid::WriteGridFaceTopology( DataBook * databook )
{
    IntField numFaceNode( this->nFace );

    for ( int iFace = 0; iFace < this->nFace; ++ iFace )
    {
        numFaceNode[ iFace ] = this->faceTopo->f2n[ iFace ].size();
    }

    ONEFLOW::HXWrite( databook, numFaceNode );

    IntField faceNodeMem;

    for ( int iFace = 0; iFace < this->nFace; ++ iFace )
    {
        int nNode = numFaceNode[ iFace ];
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            faceNodeMem.push_back( this->faceTopo->f2n[ iFace ][ iNode ] );
        }
    }
    ONEFLOW::HXWrite( databook, faceNodeMem );

    ONEFLOW::HXWrite( databook, this->faceTopo->lCell );
    ONEFLOW::HXWrite( databook, this->faceTopo->rCell );
}

void UnsGrid::WriteBoundaryTopology( DataBook * databook )
{
    int nBFace = this->faceTopo->GetNBFace();
    ONEFLOW::HXWrite( databook, nBFace );

    ONEFLOW::HXWrite( databook, this->faceTopo->bcManager->bcRecord->bcType );
    ONEFLOW::HXWrite( databook, this->faceTopo->bcManager->bcRecord->bcRegion );

    ONEFLOW::HXWrite( databook, this->interFace->nIFace );
    if ( this->interFace->nIFace > 0 )
    {
        ONEFLOW::HXWrite( databook, this->interFace->zoneId            );
        ONEFLOW::HXWrite( databook, this->interFace->localInterfaceId  );
        ONEFLOW::HXWrite( databook, this->interFace->i2b               );
    }
}

void UnsGrid::ModifyBcType( int bcType1, int bcType2 )
{
    int nBFace = this->faceTopo->bcManager->bcRecord->bcType.size();
    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bctype = this->faceTopo->bcManager->bcRecord->bcType[ iFace ];
        if ( bctype == bcType1 )
        {
            this->faceTopo->bcManager->bcRecord->bcType[ iFace ] = bcType2;
        }
    }
}

void UnsGrid::GenerateLgMapping( IFaceLink * iFaceLink )
{
    cout << "zoneIndex = " << this->id << endl;

    BcRecord * bcRecord = this->faceTopo->bcManager->bcRecord;

    this->faceTopo->bcManager->PreProcess();

    int nIFace = bcRecord->CalcNIFace();

    cout << "nIFace = " << nIFace << endl;
    this->nIFace = nIFace;
    this->interFace->Set( nIFace, this );

    if ( nIFace == 0 ) return;

    iFaceLink->Init( this );

    this->faceTopo->bcManager->bcRecord->CreateI2B( this->interFace );

    int nBFace = bcRecord->GetNBFace();

    RealField xList, yList, zList;
    IntField gINode;

    int lCount  = 0;

    for ( int iBFace = 0; iBFace < nBFace; ++ iBFace )
    {
        if ( ! BC::IsInterfaceBc( bcRecord->bcType[ iBFace ] ) )
        {
            continue;
        }
        IntField & faceNode = this->faceTopo->f2n[ iBFace ];
        int nNode = faceNode.size();

        gINode.resize( nNode );
        xList.resize( nNode );
        yList.resize( nNode );
        zList.resize( nNode );

        ONEFLOW::GetFaceCoorList( faceNode, xList, yList, zList, this->nodeMesh );
        ONEFLOW::GetCoorIdList( iFaceLink, xList, yList, zList, nNode, gINode );
        iFaceLink->CreateLink( gINode, this->id, lCount );

        ++ lCount;
    }

    cout << "local interface count = " << lCount << endl;
}

void UnsGrid::ReGenerateLgMapping( IFaceLink * iFaceLink )
{
    cout << "zoneIndex = " << this->id << endl;

    if ( ! this->faceTopo->bcManager->ExistInterface() )
    {
        return;
    }

    //modify the number of boundary faces
    //modify the number of inter faces
    //modify the face node indexes
    //modify the face node number

    this->faceTopo->ModifyFaceNodeId( iFaceLink );
    this->faceTopo->ModifyBoundaryInformation( iFaceLink );
}

void UnsGrid::UpdateOtherTopologyTerm( IFaceLink * iFaceLink )
{
    if ( ! IsValid( this->interFace ) ) return;

    this->faceTopo->UpdateOtherTopologyTerm();

    int nIFace = iFaceLink->l2g[ this->id ].size();

    this->interFace->Resize( nIFace );
    this->faceTopo->GenerateI2B( this->interFace );
}

void UnsGrid::GetMinMaxDistance( Real & dismin, Real & dismax )
{
    dismin =   LARGE;
    dismax = - LARGE;

    RealField & x = this->nodeMesh->xN;
    RealField & y = this->nodeMesh->yN;
    RealField & z = this->nodeMesh->zN;

    int nFace = this->faceTopo->GetNFace();

    Real ptTol = Tolerence::GetTol();

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        IntField & faceNode = this->faceTopo->f2n[ iFace ];
        int nNode = faceNode.size();
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            int p1 = faceNode[ iNode ];
            int p2 = faceNode[ ( iNode + 1 ) % nNode ];
            Real dx = x[ p2 ] - x[ p1 ];
            Real dy = y[ p2 ] - y[ p1 ];
            Real dz = z[ p2 ] - z[ p1 ];
            Real ds = ONEFLOW::DIST( dx, dy, dz );

            if ( ds <= ptTol ) continue;
            dismin = ONEFLOW::MIN( dismin, ds );
            dismax = ONEFLOW::MAX( dismax, ds );
        }
    }
}

void UnsGrid::WriteGrid( fstream & file )
{
    VirtualFile * vf = new VirtualFile( & file );

    vf->BeginWriteWork();

    cout << "Output nNode = " << this->nNode << " nFace = " <<  this->nFace << " nCell = " << nCell << "\n";

    ONEFLOW::HXWrite( vf, this->nNode );
    ONEFLOW::HXWrite( vf, this->nFace );
    ONEFLOW::HXWrite( vf, this->nCell );

    cout << "Output grid\n";
    ONEFLOW::HXWrite( vf, this->nodeMesh->xN );
    ONEFLOW::HXWrite( vf, this->nodeMesh->yN );
    ONEFLOW::HXWrite( vf, this->nodeMesh->zN );

    cout << "Output volBcType = " << this->volBcType << "\n";
    ONEFLOW::HXWrite( vf, this->volBcType  );

    this->WriteGridFaceTopology( vf );
    this->WriteBoundaryTopology( vf );

    vf->EndWriteWork();

    delete vf;
}

void UnsGrid::WriteGridFaceTopology( VirtualFile * vf )
{
    int nFace = this->faceTopo->f2n.size();
    IntField nFNode( nFace );

    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        nFNode[ iFace ] = this->faceTopo->f2n[ iFace ].size();
    }

    IntField fNode;
    for ( int iFace = 0; iFace < nFace; ++ iFace )
    {
        int nNode = this->faceTopo->f2n[ iFace ].size();
        for ( int iNode = 0; iNode < nNode; ++ iNode )
        {
            int nodeId = this->faceTopo->f2n[ iFace ][ iNode ];
            fNode.push_back( nodeId );
        }
    }

    ONEFLOW::HXWrite( vf, nFNode );
    ONEFLOW::HXWrite( vf, fNode );

    ONEFLOW::HXWrite( vf, this->faceTopo->lCell );
    ONEFLOW::HXWrite( vf, this->faceTopo->rCell );
}

void UnsGrid::CalcMetrics()
{
    this->AllocMetrics();

    if ( IsOneD() )
    {
        this->CalcMetrics1D();
    }
    else if ( IsTwoD() )
    {
        this->CalcMetrics2D();
    }
    else if ( IsThreeD() )
    {
        this->CalcMetrics3D();
    }
}

void UnsGrid::AllocMetrics()
{
    this->faceMesh->AllocateMetrics();
    this->cellMesh->AllocateMetrics( this->faceMesh );

}

void UnsGrid::CalcMetrics1D()
{
    this->CalcFaceCenter1D();
    this->CalcCellCenterVol1D();
    this->CalcFaceNormal1D();
    this->CalcGhostCellCenterVol1D();
}

void UnsGrid::CalcMetrics2D()
{
    this->CalcFaceNormal2D();
    this->CalcFaceCenter2D();
    this->CalcCellCenterVol2D();
}

void UnsGrid::CalcMetrics3D()
{
    this->CalcFaceNormal3D();
    this->CalcFaceCenter3D();
    this->CalcCellCenterVol3D();
}

void UnsGrid::CalcFaceCenter1D()
{
    this->faceMesh->CalcFaceCenter1D( this->nodeMesh );
}

void UnsGrid::CalcFaceNormal1D()
{
    this->faceMesh->CalcFaceNormal1D( this->nodeMesh, this->cellMesh );
}

void UnsGrid::CalcCellCenterVol1D()
{
    UInt nFace = this->faceMesh->GetNFace();
    UInt nBFace = this->faceMesh->GetNBFace();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    RealField & xcc = this->cellMesh->xcc ;
    RealField & ycc = this->cellMesh->ycc ;
    RealField & zcc = this->cellMesh->zcc ;
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

void UnsGrid::CalcGhostCellCenterVol1D()
{
    UInt nFace = this->faceMesh->GetNFace();
    UInt nBFace = this->faceMesh->GetNBFace();
    UInt numberOfCells = this->cellMesh->GetNumberOfCells();

    RealField & xcc = this->cellMesh->xcc ;
    RealField & ycc = this->cellMesh->ycc ;
    RealField & zcc = this->cellMesh->zcc ;
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

void UnsGrid::CalcFaceNormal2D()
{
    this->faceMesh->CalcFaceNormal2D( this->nodeMesh );
}

void UnsGrid::CalcFaceCenter2D()
{
    this->faceMesh->CalcFaceCenter2D( this->nodeMesh );
}

void UnsGrid::CalcCellCenterVol2D()
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

    xcc = 0;
    ycc = 0;
    zcc = 0;
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

void UnsGrid::CalcCellCenterVol3D()
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

            if ( lc == 1174 || rc == 1174 )
            {
                int kkk = 1;
            }

            xcc[ lc ] += trix;
            ycc[ lc ] += triy;
            zcc[ lc ] += triz;
            vol[ lc ] += tmp;

            xcc[ rc ] -= trix;
            ycc[ rc ] -= triy;
            zcc[ rc ] -= triz;
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

void UnsGrid::CalcFaceNormal3D()
{
    this->faceMesh->CalcFaceNormal3D( this->nodeMesh );
}

void UnsGrid::CalcFaceCenter3D()
{
    this->faceMesh->CalcFaceCenter3D( this->nodeMesh );
}

EndNameSpace