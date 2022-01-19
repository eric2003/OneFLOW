/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

#include "FaceMesh.h"
#include "NodeMesh.h"
#include "CellMesh.h"
#include "FaceTopo.h"
#include "HXMath.h"
#include "BcRecord.h"
#include <iostream>
#include <iomanip>


BeginNameSpace( ONEFLOW )

FaceMesh::FaceMesh()
{
}

FaceMesh::~FaceMesh()
{
}

HXSize_t FaceMesh::GetNFace()
{
    return faceTopo->GetNFaces();  
}

HXSize_t FaceMesh::CalcTotalFaceNodes()
{
    return faceTopo->CalcTotalFaceNodes();
}

HXSize_t FaceMesh::GetNBFace()
{
    return faceTopo->GetNBFaces();
}

void FaceMesh::SetNBFace( HXSize_t nBFaces )
{
    faceTopo->SetNBFaces( nBFaces );
}

void FaceMesh::CalcFaceCenter1D( NodeMesh * nodeMesh )
{
    HXSize_t nFaces = this->GetNFace();
    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & nodeIndex = faceTopo->faces[ iFace ];
        int p1 = nodeIndex[ 0 ];
        int p2 = nodeIndex[ 0 ];
        xfc[ iFace ] = half * ( xN[ p1 ] + xN[ p2 ] );
        yfc[ iFace ] = half * ( yN[ p1 ] + yN[ p2 ] );
        zfc[ iFace ] = half * ( zN[ p1 ] + zN[ p2 ] );
    }
}

void FaceMesh::CalcFaceNormal1D( NodeMesh * nodeMesh, CellMesh * cellMesh )
{
    HXSize_t nFaces = this->GetNFace();
    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    RealField & xcc  = cellMesh->xcc ;
    RealField & ycc  = cellMesh->ycc ;
    RealField & zcc  = cellMesh->zcc ;
    RealField & vol = cellMesh->vol;

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        int lc  = faceTopo->lCells[ iFace ];

        Real dx = xfc[ iFace ] - xcc[ lc ];
        Real dy = yfc[ iFace ] - ycc[ lc ];
        Real dz = zfc[ iFace ] - zcc[ lc ];
        Real ds = ONEFLOW::DIST( dx, dy, dz );
        
        Real factor   = 1.0 / ( ds + SMALL );
        xfn[ iFace ] = factor * dx;
        yfn[ iFace ] = factor * dy;
        zfn[ iFace ] = factor * dz;

        area[ iFace ] = 1.0;
    }
}


void FaceMesh::CalcFaceNormal2D( NodeMesh * nodeMesh )
{
    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    HXSize_t nFaces = this->GetNFace();

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & faceIndex = faceTopo->faces[ iFace ];
        int p1 = faceIndex[ 0 ];
        int p2 = faceIndex[ 1 ];

        xfn[ iFace ]  = yN[ p2 ] - yN[ p1 ];
        yfn[ iFace ]  = xN[ p1 ] - xN[ p2 ];
        zfn[ iFace ]  = 0.0;

        area[ iFace ] = ONEFLOW::DIST( xfn[ iFace ], yfn[ iFace ], zfn[ iFace ] );
    }

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        Real oArea = 1.0 / ( area[ iFace ] + SMALL );
        xfn[ iFace ] *= oArea;
        yfn[ iFace ] *= oArea;
        zfn[ iFace ] *= oArea;
    }
}

void FaceMesh::CalcFaceCenter2D( NodeMesh * nodeMesh )
{
    HXSize_t nFaces = this->GetNFace();
    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & nodeIndex = faceTopo->faces[ iFace ];
        int p1 = nodeIndex[ 0 ];
        int p2 = nodeIndex[ 1 ];
        xfc[ iFace ] = half * ( xN[ p1 ] + xN[ p2 ] );
        yfc[ iFace ] = half * ( yN[ p1 ] + yN[ p2 ] );
        zfc[ iFace ] = half * ( zN[ p1 ] + zN[ p2 ] );
        int kkk = 1;
    }
    int kkk = 1;
}

void FaceMesh::CalcFaceNormal3D( NodeMesh * nodeMesh )
{
    xfn = 0;
    yfn = 0;
    zfn = 0;

    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    HXSize_t nFaces = this->GetNFace();

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        IntField & faceIndex = faceTopo->faces[ iFace ];

        HXSize_t faceNodeNumber = faceIndex.size();
        for ( HXSize_t iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace )
        {
            int index1 = iNodeInFace;
            int index2 = ( iNodeInFace + 1 ) % faceNodeNumber;
            int p1 = faceIndex[ index1 ];
            int p2 = faceIndex[ index2 ];

            Real dx1 = xN[ p1 ];
            Real dy1 = yN[ p1 ];
            Real dz1 = zN[ p1 ];

            Real dx2 = xN[ p2 ];
            Real dy2 = yN[ p2 ];
            Real dz2 = zN[ p2 ];

            xfn[ iFace ] += half * ( dy1 * dz2 - dy2 * dz1 );
            yfn[ iFace ] += half * ( dz1 * dx2 - dz2 * dx1 );
            zfn[ iFace ] += half * ( dx1 * dy2 - dx2 * dy1 );
        }
        area[ iFace ] = ONEFLOW::DIST( xfn[ iFace ], yfn[ iFace ], zfn[ iFace ] );
    }

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        Real oArea = 1.0 / ( area[ iFace ] + SMALL );
        xfn[ iFace ] *= oArea;
        yfn[ iFace ] *= oArea;
        zfn[ iFace ] *= oArea;
    }
}

void FaceMesh::CalcFaceCenter3D( NodeMesh * nodeMesh )
{
    RealField & xN = nodeMesh->xN;
    RealField & yN = nodeMesh->yN;
    RealField & zN = nodeMesh->zN;

    HXSize_t nFaces = this->GetNFace();

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        Real x0 = 0.0;
        Real y0 = 0.0;
        Real z0 = 0.0;

        IntField & faceIndex = faceTopo->faces[ iFace ];

        HXSize_t faceNodeNumber = faceIndex.size();
        for ( HXSize_t iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace )
        {
            int index = faceIndex[ iNodeInFace ];
            x0 += xN[ index ];
            y0 += yN[ index ];
            z0 += zN[ index ];
        }

        Real factor = 1.0 / faceNodeNumber;

        x0 *= factor;
        y0 *= factor;
        z0 *= factor;

        xfc[ iFace ] = 0.0;
        yfc[ iFace ] = 0.0;
        zfc[ iFace ] = 0.0;
        Real sarea  = 0.0;

        for ( HXSize_t iNodeInFace = 0; iNodeInFace < faceNodeNumber; ++ iNodeInFace )
        {
            int index1 = iNodeInFace;
            int index2 = ( iNodeInFace + 1 ) % faceNodeNumber;
            int p1 = faceIndex[ index1 ];
            int p2 = faceIndex[ index2 ];

            Real dx1 = xN[ p1 ] - x0;
            Real dy1 = yN[ p1 ] - y0;
            Real dz1 = zN[ p1 ] - z0;

            Real dx2 = xN[ p2 ] - x0;
            Real dy2 = yN[ p2 ] - y0;
            Real dz2 = zN[ p2 ] - z0;

            Real x00 = third * ( xN[ p1 ] + xN[ p2 ] + x0 );
            Real y00 = third * ( yN[ p1 ] + yN[ p2 ] + y0 );
            Real z00 = third * ( zN[ p1 ] + zN[ p2 ] + z0 );

            Real anx = dy1 * dz2 - dy2 * dz1;
            Real any = dz1 * dx2 - dz2 * dx1;
            Real anz = dx1 * dy2 - dx2 * dy1;

            Real faceArea = half * ONEFLOW::DIST( anx, any, anz );
            Real norm = anx * xfn[ iFace ] + any * yfn[ iFace ] + anz * zfn[ iFace ];

            if ( norm < 0.0 ) faceArea  = - faceArea;
            sarea += faceArea;

            xfc[ iFace ] += faceArea * x00;
            yfc[ iFace ] += faceArea * y00;
            zfc[ iFace ] += faceArea * z00;
        }

        Real osarea = 1.0 / ( sarea + SMALL );

        xfc[ iFace ] *= osarea;
        yfc[ iFace ] *= osarea;
        zfc[ iFace ] *= osarea;
    }
}

void FaceMesh::AllocateMetrics()
{
    HXSize_t nFaces = this->GetNFace();
    this->xfc.resize( nFaces );
    this->yfc.resize( nFaces );
    this->zfc.resize( nFaces );
    this->xfn.resize( nFaces );
    this->yfn.resize( nFaces );
    this->zfn.resize( nFaces );
    this->area.resize( nFaces );
    this->vfx.resize( nFaces );
    this->vfy.resize( nFaces );
    this->vfz.resize( nFaces );
    this->vfn.resize( nFaces );
    this->vfx = 0;
    this->vfy = 0;
    this->vfz = 0;
    this->vfn = 0;
    HXSize_t nBFaces = this->GetNBFace();
    this->faceTopo->bcManager->bcRecord->bcType.resize( nBFaces );
}

EndNameSpace
