/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

#include "FaceSolver.h"
#include "ElementHome.h"
#include "Stop.h"
#include "FaceTopo.h"
#include "CgnsSection.h"
#include <iostream>
#include <algorithm>


BeginNameSpace( ONEFLOW )

FaceSolver::FaceSolver()
{
    this->faceBcKey = new IntField();
    this->faceBcType = new IntField();
    this->childFid = new LinkField();

    this->refFaces = new std::set< HXMid<int> >();
    this->faceTopo = new FaceTopo();
}

FaceSolver::~FaceSolver()
{
    delete this->faceBcKey;
    delete this->faceBcType;
    delete this->childFid;

    delete this->refFaces;
    delete this->faceTopo;
}

int FaceSolver::FindFace( HXMid<int> & face )
{
    std::set< HXMid<int> >::iterator iter = this->refFaces->find( face );
    if ( iter == this->refFaces->end() )
    {
        return ONEFLOW::INVALID_INDEX;
    }
    return iter->id;
}

bool FaceSolver::CheckBcFace( IntSet & bcVertex, IntField & nodeId )
{
    int size = nodeId.size();
    for ( int iNode = 0; iNode < size; ++ iNode )
    {
        IntSet::iterator iter = bcVertex.find( nodeId[ iNode ] );
        if ( iter == bcVertex.end() )
        {
            return false;
        }
    }
    return true;
}

void FaceSolver::ScanPolygonFace( CgnsSection * cgnsSection )
{
    std::vector<int> faceNodes;
    for ( int iElem = 0; iElem < cgnsSection->nElement; ++ iElem )
    {
        int st = cgnsSection->ePosList[ iElem ];
        int ed = cgnsSection->ePosList[ iElem + 1 ];
        int nNode = ed - st;
        faceNodes.resize( 0 );
        for ( int i = st; i < ed; ++ i )
        {
            int node = cgnsSection->connList[ i ];
            faceNodes.push_back( node );
        }
        HXMid<int> fMid( nNode, this->refFaces->size() );
        fMid.data = faceNodes;
        std::sort( fMid.data.begin(), fMid.data.end() );
        int gFid = this->FindFace( fMid );
        if ( gFid == ONEFLOW::INVALID_INDEX )
        {
            this->refFaces->insert( fMid );
            this->faceTopo->faces.push_back( faceNodes );
            this->faceTopo->fTypes.push_back(cgnsSection->eType);
            this->faceTopo->faceFlags.push_back(0);
        }
        int kkk = 1;
    }
    int kkk = 1;
}

void FaceSolver::ResizeAll()
{
    this->faceTopo->ResizeAll();
    int nFaces = this->faceTopo->faces.size();
    this->faceBcType->resize( nFaces );
    this->faceBcKey->resize( nFaces );
    this->childFid->resize( nFaces );
}

void FaceSolver::ScanPolyhedronElement( CgnsSection * cgnsSection )
{
    std::vector<int> faceIds;
    for ( int iElem = 0; iElem < cgnsSection->nElement; ++ iElem )
    {
        int st = cgnsSection->ePosList[ iElem ];
        int ed = cgnsSection->ePosList[ iElem + 1 ];
        int nFace = ed - st;
        faceIds.resize( 0 );
        for ( int i = st; i < ed; ++ i )
        {
            int polygonFaceId = std::abs(cgnsSection->connList[ i ]);
            faceIds.push_back( polygonFaceId );

            int faceFlags = this->faceTopo->faceFlags[ polygonFaceId ];

            if ( faceFlags == 0 ) //face left element not set
            {
                this->ResizeAll();
                this->faceTopo->faceFlags[ polygonFaceId ] = 1;
                this->faceTopo->lCells[ polygonFaceId ] = iElem;
                this->faceTopo->rCells[ polygonFaceId ] = ONEFLOW::INVALID_INDEX;

                (*this->faceBcType)[ polygonFaceId ] = ONEFLOW::INVALID_INDEX;
                (*this->faceBcKey )[ polygonFaceId ] = ONEFLOW::INVALID_INDEX;
            }
            else
            {
                this->faceTopo->rCells[ polygonFaceId ] = iElem;
            }

        }
        int kkk = 1;
    }
    int kkk = 1;
}

void FaceSolver::ScanElementFace( CgIntField & eNodeId, int eType, int eId )
{
    UnitElement * unitElement = ElementHome::GetUnitElement( eType );

    //composite Element not to be involved in analysis !!!
    int nElemFace = unitElement->faceList.size();
    for ( int iFace = 0; iFace < nElemFace; ++ iFace )
    {
        IntField & rNodeId = unitElement->faceList[ iFace ];
        int fType = unitElement->GetFaceType( iFace );
         
        int nNodes = rNodeId.size();

        IntField aNodeId;
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            aNodeId.push_back( eNodeId[ rNodeId[ iNode ] ] );
        }                                                              

        HXMid<int> fMid( nNodes, this->faceTopo->faces.size() );
        fMid.data = aNodeId;
        std::sort( fMid.data.begin(), fMid.data.end() );
        int gFid = this->FindFace( fMid );

        if ( gFid == ONEFLOW::INVALID_INDEX )
        {
            int totalfn    = this->refFaces->size();
            int faceNumber = this->faceTopo->lCells.size();
            if ( totalfn != faceNumber )
            {
                std::cout << "totalfn != faceNumber " << totalfn << " " << faceNumber << std::endl;
                Stop("");
            }
            this->refFaces->insert( fMid );

            this->faceTopo->lCells.push_back( eId );
            this->faceTopo->rCells.push_back( ONEFLOW::INVALID_INDEX );

            this->faceBcType->push_back( ONEFLOW::INVALID_INDEX );
            this->faceBcKey->push_back( ONEFLOW::INVALID_INDEX );
            this->faceTopo->fTypes.push_back( fType );

            this->faceTopo->faces.push_back( aNodeId );
            this->childFid->resize( totalfn + 1 );
        }
        else
        {
            if ((this->faceTopo->lCells)[gFid] == ONEFLOW::INVALID_INDEX)
            {
                //This shows that although this aspect exists, it has not been dealt with due to various reasons
                (this->faceTopo->lCells)[gFid] = eId; //For example, a new volume element surface is added during the splitting process
            }
            else
            {
                if ( (this->faceTopo->rCells)[gFid] == ONEFLOW::INVALID_INDEX )
                {
                    if ((this->faceTopo->lCells)[gFid] != eId)
                    {
                        (this->faceTopo->rCells)[gFid] = eId;
                    }
                }
            }
        }
    }
}

void FaceSolver::ScanBcFace( IntSet& bcVertex, int bcType, int bcNameId )
{
    int nBFaces = 0;

    //std::cout << " this->faceTopo = " << this->faceTopo << "\n";
    int nFaces = this->faceTopo->lCells.size();

    std::cout << " nFaces = " << nFaces << "\n";
    int nTraditionalBc = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rCell = ( this->faceTopo->rCells )[ iFace ];

        if ( rCell == ONEFLOW::INVALID_INDEX )
        {
            ++ nTraditionalBc;
        }
    }
    std::cout << " nTraditionalBc = " << nTraditionalBc << "\n";


    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( iFace % 200000 == 0 ) 
        {
            //std::cout << " iFace = " << iFace << " numberOfTotalFaces = " << nFaces << std::endl;
        }
        int originalBcType = ( * this->faceBcType )[ iFace ];
        int rCell     = ( this->faceTopo->rCells )[ iFace ];

        if ( ( rCell          == ONEFLOW::INVALID_INDEX ) && 
             ( originalBcType == ONEFLOW::INVALID_INDEX ) )
        {
            if ( this->CheckBcFace( bcVertex, ( this->faceTopo->faces )[ iFace ] ) )
            {
                ++ nBFaces;

                ( * this->faceBcType )[ iFace ] = bcType;
                ( * this->faceBcKey  )[ iFace ] = bcNameId;
            }
        }
    }

    //std::cout << " nBFaces = " << nBFaces << std::endl;
    int kkk = 1;
}

void FaceSolver::ScanBcFaceDetail( IntSet& bcVertex, int bcType, int bcNameId )
{
    int nFaces = this->faceTopo->lCells.size();
    std::cout << " nFaces = " << nFaces << "\n";

    int nTraditionalBc = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rCell = ( this->faceTopo->rCells )[ iFace ];

        if ( rCell == ONEFLOW::INVALID_INDEX )
        {
            ++ nTraditionalBc;
        }
    }
    std::cout << " nTraditionalBc = " << nTraditionalBc << "\n";

    int nBFaces = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( iFace % 200000 == 0 ) 
        {
            //std::cout << " iFace = " << iFace << " numberOfTotalFaces = " << nFaces << std::endl;
        }
        int originalBcType = ( * this->faceBcType )[ iFace ];

        if ( originalBcType == ONEFLOW::INVALID_INDEX )
        {
            if ( this->CheckBcFace( bcVertex, ( this->faceTopo->faces )[ iFace ] ) )
            {
                ++ nBFaces;

                ( * this->faceBcType )[ iFace ] = bcType;
                ( * this->faceBcKey  )[ iFace ] = bcNameId;
            }
        }
    }

    std::cout << " nFinalBcFace = " << nBFaces << " bcType = " << bcType << std::endl;
    int kkk = 1;
}

void FaceSolver::ScanInterfaceBc()
{
    int nFaces = this->faceTopo->lCells.size();

    int bcNameId = -1;
    int nInterFace = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int originalBcType = ( * this->faceBcType )[ iFace ];

        if ( originalBcType == ONEFLOW::INVALID_INDEX )
        {
            nInterFace ++;
            ( * this->faceBcType )[ iFace ] = BCTypeNull;
            ( * this->faceBcKey  )[ iFace ] = bcNameId;
        }
    }

    std::cout << " nInterFace = " << nInterFace << std::endl;
}

int FaceSolver::GetNSimpleFace()
{
    int nSimpleFace = 0;

    for ( int iFace = 0; iFace < this->faceTopo->faces.size(); ++ iFace )
    {
        int nCFace = ( * this->childFid )[ iFace ].size();
        if ( nCFace == 0 )
        {
            ++ nSimpleFace;
        }
    }
    return nSimpleFace;
}


EndNameSpace
