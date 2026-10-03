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

#include "UnsGrid.h"
#include <memory>
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
#include "IFaceLink.h"
#include "Boundary.h"
#include "Tolerance.h"
#include "Dimension.h"
#include "LogFile.h"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <iomanip>


BeginNameSpace( ONEFLOW )

REGISTER_GRID( UnsGrid )

UnsGrid * UnsGridCast( Grid * gridIn )
{
    return static_cast< UnsGrid * >( gridIn );
}

UnsGrid::UnsGrid() = default;

UnsGrid::~UnsGrid() = default;

void UnsGrid::Init()
{
    this->cellMesh.reset();
    this->faceMesh.reset();
    this->faceTopo.reset();

    this->BasicInit();
    this->faceMesh = std::make_unique< FaceMesh >();
    this->cellMesh = std::make_unique< CellMesh >();
    this->SetFaceTopo( std::make_unique< FaceTopo >() );
}

FaceTopo & UnsGrid::GetFaceTopo()
{
    return *this->faceTopo;
}

const FaceTopo & UnsGrid::GetFaceTopo() const
{
    return *this->faceTopo;
}

FaceMesh & UnsGrid::GetFaceMesh()
{
    return *this->faceMesh;
}

const FaceMesh & UnsGrid::GetFaceMesh() const
{
    return *this->faceMesh;
}

CellMesh & UnsGrid::GetCellMesh()
{
    return *this->cellMesh;
}

const CellMesh & UnsGrid::GetCellMesh() const
{
    return *this->cellMesh;
}

void UnsGrid::SetFaceTopo( std::unique_ptr< FaceTopo > faceTopo )
{
    this->faceTopo = std::move( faceTopo );
    this->GetFaceTopo().BindGrid( *this );
    this->GetFaceMesh().BindFaceTopo( this->GetFaceTopo() );
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
    std::cout << "Reading unstructured grid data files......\n";
    //Read the number of nodes, number of elements surface and number of elements

    ONEFLOW::HXRead( databook, this->nNodes );
    ONEFLOW::HXRead( databook, this->nFaces );
    ONEFLOW::HXRead( databook, this->nCells );

    std::cout << "Grid dimension = " << Dim::dimension << std::endl;

    std::cout << " number of nodes    : " << this->nNodes << std::endl;
    std::cout << " number of surfaces : " << this->nFaces << std::endl;
    std::cout << " number of elements : " << this->nCells << std::endl;

    this->nodeMesh->CreateNodes( this->nNodes );
    this->GetCellMesh().GetCellTopo().Alloc( this->nCells );

    ONEFLOW::HXRead( databook, this->nodeMesh->xN );
    ONEFLOW::HXRead( databook, this->nodeMesh->yN );
    ONEFLOW::HXRead( databook, this->nodeMesh->zN );

    std::cout << "The grid nodes have been read\n";
    ONEFLOW::HXRead( databook, this->volBcType  );

    this->nodeMesh->CalcMinMaxBox();
    this->ReadGridFaceTopology( databook );
    this->ReadBoundaryTopology( databook );
    this->NormalizeBc();

    std::cout << "All the computing information is ready!\n";
}

void UnsGrid::NormalizeBc()
{
    for ( int iFace = 0; iFace < this->nBFaces; ++ iFace )
    {
        this->GetFaceTopo().GetRightCells()[ iFace ] = iFace + this->nCells;
    }
}

void UnsGrid::ReadGridFaceTopology( DataBook * databook )
{
    this->GetFaceTopo().GetFaces().resize( this->nFaces );
    this->GetFaceTopo().GetLeftCells().resize( this->nFaces );
    this->GetFaceTopo().GetRightCells().resize( this->nFaces );
    this->GetFaceTopo().GetFaceTypes().resize( this->nFaces );

    IntField numFaceNode( this->nFaces );

    ONEFLOW::HXRead( databook, numFaceNode );

    int nsum = ONEFLOW::SUM( numFaceNode );

    std::cout << "Setting the connection mode of face to point......\n";
    IntField faceNodeMem( nsum );

    ONEFLOW::HXRead( databook, faceNodeMem );

    int ipos = 0;
    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        int nNodes = numFaceNode[ iFace ];
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            int pid = faceNodeMem[ ipos ++ ];
            this->GetFaceTopo().GetFaces()[ iFace ].push_back( pid );
        }
    }

    std::cout << "Setting the connection mode of face to cell......\n";

    ONEFLOW::HXRead( databook, this->GetFaceTopo().GetLeftCells() );
    ONEFLOW::HXRead( databook, this->GetFaceTopo().GetRightCells() );

    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        if ( this->GetFaceTopo().GetLeftCells()[ iFace ] < 0 )
        {
            //need to reverse the node ordering
            IntField & f2n = this->GetFaceTopo().GetFaces()[ iFace ];
            std::reverse( f2n.begin(), f2n.end() );
            // now reverse leftCellIndex  and rightCellIndex
            ONEFLOW::SWAP( this->GetFaceTopo().GetLeftCells()[ iFace ], this->GetFaceTopo().GetRightCells()[ iFace ] );
        }
    }
}

void UnsGrid::ReadBoundaryTopology( DataBook * databook )
{
    std::cout << "Setting the boundary condition......\n";
    ONEFLOW::HXRead( databook, this->nBFaces );
    this->GetFaceTopo().SetNBFaces( this->nBFaces );

    //std::cout << " nBFaces = " << this->nBFaces << std::endl;

    //Setting boundary conditions
    BcRecord * bcRecord = &this->GetFaceTopo().GetBcRecord();
    ONEFLOW::HXRead( databook, bcRecord->bcType );
    ONEFLOW::HXRead( databook, bcRecord->bcNameId );
    ONEFLOW::HXRead( databook, this->nIFaces );
    std::cout << " nBFaces = " << this->nBFaces;
    std::cout << " nIFaces = " << this->nIFaces << std::endl;
    this->interFace->Set( this->nIFaces, this );

    if ( this->nIFaces > 0 )
    {
        ONEFLOW::HXRead( databook, this->interFace->zoneId            );
        ONEFLOW::HXRead( databook, this->interFace->localInterfaceId  );
        ONEFLOW::HXRead( databook, this->interFace->i2b               );
    }
}

void UnsGrid::WriteGrid( DataBook * databook )
{
    std::cout << "Grid dimension = " << Dim::dimension << std::endl;

    std::cout << " number of nodes    : " << this->nNodes << std::endl;
    std::cout << " number of surfaces : " << this->nFaces << std::endl;
    std::cout << " number of elements : " << this->nCells << std::endl;

    ONEFLOW::HXWrite( databook, this->nNodes );
    ONEFLOW::HXWrite( databook, this->nFaces );
    ONEFLOW::HXWrite( databook, this->nCells );

    ONEFLOW::HXWrite( databook, this->nodeMesh->xN );
    ONEFLOW::HXWrite( databook, this->nodeMesh->yN );
    ONEFLOW::HXWrite( databook, this->nodeMesh->zN );

    ONEFLOW::HXWrite( databook, this->volBcType  );

    if ( ONEFLOW::IsOneD() )
    {
        this->WriteGridFaceTopology1D( databook );
        this->WriteBoundaryTopology1D( databook );
    }
    else
    {
        this->WriteGridFaceTopology( databook );
        this->WriteBoundaryTopology( databook );
    }
}

void UnsGrid::WriteGridFaceTopology1D( DataBook * databook )
{
    std::cout << " Reading eTypes\n";

    //write element types
    int ntmpElements = this->GetCellMesh().GetCellTopo().eTypes.size();
    ONEFLOW::HXWrite( databook, this->GetCellMesh().GetCellTopo().eTypes );

    //write face types
    int ntmpFaces = this->GetFaceTopo().GetFaceTypes().size();
    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetFaceTypes() );

    IntField numFaceNode( this->nFaces );

    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        numFaceNode[ iFace ] = this->GetFaceTopo().GetFaces()[ iFace ].size();
    }

    ONEFLOW::HXWrite( databook, numFaceNode );

    IntField faceNodeMem;

    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        int nNodes = numFaceNode[ iFace ];
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            faceNodeMem.push_back( this->GetFaceTopo().GetFaces()[ iFace ][ iNode ] );
        }
    }
    ONEFLOW::HXWrite( databook, faceNodeMem );

    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetLeftCells() );
    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetRightCells() );
}

void UnsGrid::WriteGridFaceTopology( DataBook * databook )
{
    IntField numFaceNode( this->nFaces );

    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        numFaceNode[ iFace ] = this->GetFaceTopo().GetFaces()[ iFace ].size();
    }

    ONEFLOW::HXWrite( databook, numFaceNode );

    IntField faceNodeMem;

    for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
    {
        int nNodes = numFaceNode[ iFace ];
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            faceNodeMem.push_back( this->GetFaceTopo().GetFaces()[ iFace ][ iNode ] );
        }
    }
    ONEFLOW::HXWrite( databook, faceNodeMem );

    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetLeftCells() );
    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetRightCells() );
}

void UnsGrid::WriteBoundaryTopology( DataBook * databook )
{
    int nBFaces = this->GetFaceTopo().GetNBFaces();
    ONEFLOW::HXWrite( databook, nBFaces );

    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetBcRecord().bcType );
    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetBcRecord().bcNameId );

    ONEFLOW::HXWrite( databook, this->interFace->nIFaces );
    if ( this->interFace->nIFaces > 0 )
    {
        ONEFLOW::HXWrite( databook, this->interFace->zoneId            );
        ONEFLOW::HXWrite( databook, this->interFace->localInterfaceId  );
        ONEFLOW::HXWrite( databook, this->interFace->i2b               );
    }
}

void UnsGrid::WriteBoundaryTopology1D( DataBook * databook )
{
    int nBFaces = this->GetFaceTopo().GetNBFaces();
    ONEFLOW::HXWrite( databook, nBFaces );

    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetBcRecord().bcType );
    ONEFLOW::HXWrite( databook, this->GetFaceTopo().GetBcRecord().bcNameId );

    ONEFLOW::HXWrite( databook, this->interFace->nIFaces );
    if ( this->interFace->nIFaces > 0 )
    {
        ONEFLOW::HXWrite( databook, this->interFace->zoneId            );
        ONEFLOW::HXWrite( databook, this->interFace->localInterfaceId  );
        ONEFLOW::HXWrite( databook, this->interFace->i2b               );
    }
}


void UnsGrid::ModifyBcType( int bcType1, int bcType2 )
{
    int nBFaces = this->GetFaceTopo().GetBcRecord().bcType.size();
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bctype = this->GetFaceTopo().GetBcRecord().bcType[ iFace ];
        if ( bctype == bcType1 )
        {
            this->GetFaceTopo().GetBcRecord().bcType[ iFace ] = bcType2;
        }
    }
}

void UnsGrid::GenerateLgMapping( IFaceLink * iFaceLink )
{
    std::cout << "zoneIndex = " << this->id << std::endl;

    BcRecord * bcRecord = &this->GetFaceTopo().GetBcRecord();

    this->GetFaceTopo().PrepareBoundaryConditions();

    int nIFaces = bcRecord->CalcNIFace();

    std::cout << "nIFaces = " << nIFaces << std::endl;
    this->nIFaces = nIFaces;
    this->interFace->Set( nIFaces, this );

    if ( nIFaces == 0 ) return;

    iFaceLink->Init( this );

    this->GetFaceTopo().GenerateI2B( this->interFace.get() );

    int nBFaces = bcRecord->GetNBFace();

    RealField xList, yList, zList;
    IntField gINode;

    int lCount  = 0;

    for ( int iBFace = 0; iBFace < nBFaces; ++ iBFace )
    {
        if ( ! BC::IsInterfaceBc( bcRecord->bcType[ iBFace ] ) )
        {
            continue;
        }
        IntField & faceNode = this->GetFaceTopo().GetFaces()[ iBFace ];
        int nNodes = faceNode.size();

        gINode.resize( nNodes );
        xList.resize( nNodes );
        yList.resize( nNodes );
        zList.resize( nNodes );

        ONEFLOW::GetFaceCoorList( faceNode, xList, yList, zList, this->nodeMesh.get() );
        ONEFLOW::GetCoorIdList( iFaceLink, xList, yList, zList, nNodes, gINode );
        iFaceLink->CreateLink( gINode, this->id, lCount );

        ++ lCount;
    }

    std::cout << "local interface count = " << lCount << std::endl;
}

void UnsGrid::ReGenerateLgMapping( IFaceLink * iFaceLink )
{
    std::cout << "zoneIndex = " << this->id << std::endl;

    if ( ! this->GetFaceTopo().HasInterfaceBoundary() )
    {
        return;
    }

    //modify the number of boundary faces
    //modify the number of inter faces
    //modify the face node indexes
    //modify the face node number

    this->GetFaceTopo().ModifyFaceNodeId( iFaceLink );
    this->GetFaceTopo().ModifyBoundaryInformation( iFaceLink );
}

void UnsGrid::UpdateOtherTopologyTerm( IFaceLink * iFaceLink )
{
    if ( ! IsValid( this->interFace.get() ) ) return;

    this->GetFaceTopo().UpdateOtherTopologyTerm();

    int nIFaces = iFaceLink->l2g[ this->id ].size();

    this->interFace->Resize( nIFaces );
    this->GetFaceTopo().GenerateI2B( this->interFace.get() );
}

void UnsGrid::GetMinMaxDistance( Real & dismin, Real & dismax )
{
    dismin =   LARGE;
    dismax = - LARGE;

    RealField & x = this->nodeMesh->xN;
    RealField & y = this->nodeMesh->yN;
    RealField & z = this->nodeMesh->zN;

    int nFaces = this->GetFaceTopo().GetNFaces();

    Real ptTol = Tolerance::GetTol();

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & faceNode = this->GetFaceTopo().GetFaces()[ iFace ];
        int nNodes = faceNode.size();
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            int p1 = faceNode[ iNode ];
            int p2 = faceNode[ ( iNode + 1 ) % nNodes ];
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

void UnsGrid::WriteGrid( std::fstream & file )
{
    auto databook = std::make_unique< DataBook >();
    this->WriteGrid( databook.get() );
    databook->WriteFile( file );
}


void UnsGrid::CalcMetrics()
{
    this->AllocMetrics();

    //if ( this->IsOneD() )
    //{
    //    this->CalcMetrics1D();
    //}
    //else if ( this->IsTwoD() )
    //{
    //    this->CalcMetrics2D();
    //}
    //else if ( this->IsThreeD() )
    //{
    //    this->CalcMetrics3D();
    //}

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

void UnsGrid::AllocMetrics()
{
    this->GetFaceMesh().AllocateMetrics();
    this->GetCellMesh().AllocateMetrics( this->GetFaceMesh() );
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
    this->GetFaceMesh().CalcFaceCenter1D( this->nodeMesh.get() );
}

void UnsGrid::CalcFaceNormal1D()
{
    this->GetFaceMesh().CalcFaceNormal1D( this->nodeMesh.get(), &this->GetCellMesh() );
}

void UnsGrid::CalcCellCenterVol1D()
{
    HXSize_t nFaces = this->GetFaceMesh().GetNFace();
    HXSize_t nBFaces = this->GetFaceMesh().GetNBFace();
    HXSize_t numberOfCells = this->GetCellMesh().GetNumberOfCells();

    RealField & xcc = this->GetCellMesh().xcc ;
    RealField & ycc = this->GetCellMesh().ycc ;
    RealField & zcc = this->GetCellMesh().zcc ;
    RealField & vol = this->GetCellMesh().vol;

    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    CellTopo & cellTopo = this->GetCellMesh().GetCellTopo();
    FaceTopo & faceTopo = this->GetFaceMesh().GetFaceTopo();

    for ( HXSize_t iCell = 0; iCell < numberOfCells; ++ iCell )
    {
        IntField & element = cellTopo.elements[ iCell ];
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
    HXSize_t nFaces = this->GetFaceMesh().GetNFace();
    HXSize_t nBFaces = this->GetFaceMesh().GetNBFace();
    HXSize_t numberOfCells = this->GetCellMesh().GetNumberOfCells();

    RealField & xcc = this->GetCellMesh().xcc ;
    RealField & ycc = this->GetCellMesh().ycc ;
    RealField & zcc = this->GetCellMesh().zcc ;
    RealField & vol = this->GetCellMesh().vol;

    RealField & xfn = this->GetFaceMesh().xfn;
    RealField & yfn = this->GetFaceMesh().yfn;
    RealField & zfn = this->GetFaceMesh().zfn;

    RealField & xfc = this->GetFaceMesh().xfc;
    RealField & yfc = this->GetFaceMesh().yfc;
    RealField & zfc = this->GetFaceMesh().zfc;

    RealField & area = this->GetFaceMesh().area;

    CellTopo & cellTopo = this->GetCellMesh().GetCellTopo();
    FaceTopo & faceTopo = this->GetFaceMesh().GetFaceTopo();

    // For ghost cells
    for ( HXSize_t iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int lc  = faceTopo.GetLeftCells()[ iFace ];
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
    this->GetFaceMesh().CalcFaceNormal2D( this->nodeMesh.get() );
}

void UnsGrid::CalcFaceCenter2D()
{
    this->GetFaceMesh().CalcFaceCenter2D( this->nodeMesh.get() );
}

void UnsGrid::CalcCellCenterVol2D()
{
    HXSize_t nFaces = this->GetFaceMesh().GetNFace();
    HXSize_t nBFaces = this->GetFaceMesh().GetNBFace();
    HXSize_t numberOfCells = this->GetCellMesh().GetNumberOfCells();

    RealField & xcc  = this->GetCellMesh().xcc ;
    RealField & ycc  = this->GetCellMesh().ycc ;
    RealField & zcc  = this->GetCellMesh().zcc ;
    RealField & vol = this->GetCellMesh().vol;

    RealField & xfn = this->GetFaceMesh().xfn;
    RealField & yfn = this->GetFaceMesh().yfn;
    RealField & zfn = this->GetFaceMesh().zfn;

    RealField & xfc = this->GetFaceMesh().xfc;
    RealField & yfc = this->GetFaceMesh().yfc;
    RealField & zfc = this->GetFaceMesh().zfc;

    RealField & area = this->GetFaceMesh().area;

    CellTopo & cellTopo = this->GetCellMesh().GetCellTopo();
    FaceTopo & faceTopo = this->GetFaceMesh().GetFaceTopo();

    xcc = 0;
    ycc = 0;
    zcc = 0;
    vol = 0;

    for ( HXSize_t iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int lc = faceTopo.GetLeftCells()[ iFace ];
        Real dot = ( xfc[ iFace ] * xfn[ iFace ] +
                     yfc[ iFace ] * yfn[ iFace ] +
                     zfc[ iFace ] * zfn[ iFace ] ) * area[ iFace ];
        xcc [ lc ] += xfc[ iFace ] * dot;
        ycc [ lc ] += yfc[ iFace ] * dot;
        zcc [ lc ] += zfc[ iFace ] * dot;
        vol[ lc ] += dot;
    }

    // For interior cell faces
    for ( HXSize_t iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        int lc = faceTopo.GetLeftCells()[ iFace ];
        int rc = faceTopo.GetRightCells()[ iFace ];
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

    HXSize_t numberOfCellsHaveNegativeVolumes = 0;
    Real minvol = LARGE, maxvol = 0.0;
    HXSize_t indexMinv = 0, indexMaxv = 0;
    for ( HXSize_t iCell = 0; iCell < numberOfCells; ++ iCell )
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
    for ( HXSize_t iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int lc = faceTopo.GetLeftCells()[ iFace ];
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
    HXSize_t nFaces = this->GetFaceMesh().GetNFace();
    HXSize_t nBFaces = this->GetFaceMesh().GetNBFace();
    HXSize_t numberOfCells = this->GetCellMesh().GetNumberOfCells();

    RealField & xcc  = this->GetCellMesh().xcc ;
    RealField & ycc  = this->GetCellMesh().ycc ;
    RealField & zcc  = this->GetCellMesh().zcc ;
    RealField & vol = this->GetCellMesh().vol;

    RealField & xfn = this->GetFaceMesh().xfn;
    RealField & yfn = this->GetFaceMesh().yfn;
    RealField & zfn = this->GetFaceMesh().zfn;

    RealField & xfc = this->GetFaceMesh().xfc;
    RealField & yfc = this->GetFaceMesh().yfc;
    RealField & zfc = this->GetFaceMesh().zfc;

    RealField & area = this->GetFaceMesh().area;

    CellTopo & cellTopo = this->GetCellMesh().GetCellTopo();
    FaceTopo & faceTopo = this->GetFaceMesh().GetFaceTopo();

    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    xcc  = 0;
    ycc  = 0;
    zcc  = 0;
    vol = 0;

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        int lc = faceTopo.GetLeftCells()[ iFace ];
        int rc = faceTopo.GetRightCells()[iFace];

        IntField & faceIndex = faceTopo.GetFaces()[ iFace ];

        HXSize_t faceNodeNumber = faceIndex.size();
        for ( HXSize_t iNode = 0; iNode < faceNodeNumber; ++ iNode )
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

    HXSize_t cell = 0;
    Real minvol = LARGE, maxvol = 0.0;
    int indexMinv = 0, indexMaxv = 0;
    for ( HXSize_t iCell = 0; iCell < numberOfCells; ++ iCell )
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

    if ( cell ) std::cout << cell << " cells have negative vols \n";

    // For ghost cells
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int lc = faceTopo.GetLeftCells()[ iFace ];
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
    this->GetFaceMesh().CalcFaceNormal3D( this->nodeMesh.get() );
}

void UnsGrid::CalcFaceCenter3D()
{
    this->GetFaceMesh().CalcFaceCenter3D( this->nodeMesh.get() );
}

EndNameSpace
