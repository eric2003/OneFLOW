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

#include "FaceSolver.h"
#include "ElementHome.h"
#include "Fatal.h"
#include "FaceTopo.h"
#include "CgnsSection.h"
#include <iostream>
#include <algorithm>


BeginNameSpace( ONEFLOW )

// FaceSolver.cpp (Constructors and Destructors)
FaceSolver::FaceSolver()
{
    // [Refactored] faceBcKey, faceBcType, childFid are now value types, 
    // no need to 'new' them.
    this->faceTopo = std::make_unique< FaceTopo >();
}

FaceSolver::~FaceSolver()
    = default;

FaceTopo & FaceSolver::GetFaceTopo()
{
    return *this->faceTopo;
}

const FaceTopo & FaceSolver::GetFaceTopo() const
{
    return *this->faceTopo;
}

std::unique_ptr< FaceTopo > FaceSolver::TakeFaceTopo() noexcept
{
    return std::move( this->faceTopo );
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
    //std::vector<int> faceNodes;
    IntField faceNodes;
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

        auto [faceIndex, isNew] = faceLookup.FindOrAdd(faceNodes);
        if ( isNew )
        {
            // New face: ID is set to the current number of faces. 
            int newId = static_cast<int>(this->faceTopo->GetFaces().size());
            this->faceTopo->GetFaces().push_back(faceNodes);           // Preserve original order
            this->faceTopo->GetFaceTypes().push_back(cgnsSection->eType);
            this->faceTopo->GetFaceFlags().push_back(0);
        }
    }
}

void FaceSolver::ResizeAll()
{
    this->faceTopo->ResizeAll();
    int nFaces = this->faceTopo->GetFaces().size();
    this->faceBcType.resize( nFaces );
    this->faceBcKey.resize( nFaces );
    this->childFid.resize( nFaces );
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

            int faceFlags = this->faceTopo->GetFaceFlags()[ polygonFaceId ];

            if ( faceFlags == 0 ) //face left element not set
            {
                this->ResizeAll();
                this->faceTopo->GetFaceFlags()[ polygonFaceId ] = 1;
                this->faceTopo->GetLeftCells()[ polygonFaceId ] = iElem;
                this->faceTopo->GetRightCells()[ polygonFaceId ] = ONEFLOW::INVALID_INDEX;

                this->faceBcType[ polygonFaceId ] = ONEFLOW::INVALID_INDEX;
                this->faceBcKey[ polygonFaceId ] = ONEFLOW::INVALID_INDEX;
            }
            else
            {
                this->faceTopo->GetRightCells()[ polygonFaceId ] = iElem;
            }

        }

    }

}

void FaceSolver::ScanElementFace( CgIntField & eNodeId, int eType, int eId )
{
    UnitElement & unitElement = ElementHome::GetUnitElement( eType );

    //composite Element not to be involved in analysis !!!
    int nElemFace = unitElement.faceList.size();
    for ( int iFace = 0; iFace < nElemFace; ++ iFace )
    {
        IntField & rNodeId = unitElement.faceList[ iFace ];
        int fType = unitElement.GetFaceType( iFace );
         
        int nNodes = rNodeId.size();

        IntField aNodeId;
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            aNodeId.push_back( eNodeId[ rNodeId[ iNode ] ] );
        }                                                              

        auto [gFid, isNew]  = this->faceLookup.FindOrAdd( aNodeId );

        if ( isNew )
        {
            int totalfn = this->faceLookup.Size();

            this->faceTopo->GetLeftCells().push_back( eId );
            this->faceTopo->GetRightCells().push_back( ONEFLOW::INVALID_INDEX );

            this->faceBcType.push_back( ONEFLOW::INVALID_INDEX );
            this->faceBcKey.push_back( ONEFLOW::INVALID_INDEX );
            this->faceTopo->GetFaceTypes().push_back( fType );

            this->faceTopo->GetFaces().push_back( aNodeId );
            this->childFid.resize( totalfn );
        }
        else
        {
            if ((this->faceTopo->GetLeftCells())[gFid] == ONEFLOW::INVALID_INDEX)
            {
                //This shows that although this aspect exists, it has not been dealt with due to various reasons
                (this->faceTopo->GetLeftCells())[gFid] = eId; //For example, a new volume element surface is added during the splitting process
            }
            else
            {
                if ( (this->faceTopo->GetRightCells())[gFid] == ONEFLOW::INVALID_INDEX )
                {
                    if ((this->faceTopo->GetLeftCells())[gFid] != eId)
                    {
                        (this->faceTopo->GetRightCells())[gFid] = eId;
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
    int nFaces = this->faceTopo->GetLeftCells().size();

    std::cout << " nFaces = " << nFaces << "\n";
    int nTraditionalBc = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rCell = ( this->faceTopo->GetRightCells() )[ iFace ];

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
        int originalBcType = this->faceBcType[ iFace ];
        int rCell     = ( this->faceTopo->GetRightCells() )[ iFace ];

        if ( ( rCell          == ONEFLOW::INVALID_INDEX ) && 
             ( originalBcType == ONEFLOW::INVALID_INDEX ) )
        {
            if ( this->CheckBcFace( bcVertex, ( this->faceTopo->GetFaces() )[ iFace ] ) )
            {
                ++ nBFaces;

                this->faceBcType[ iFace ] = bcType;
                this->faceBcKey[ iFace ] = bcNameId;
            }
        }
    }

    //std::cout << " nBFaces = " << nBFaces << std::endl;

}

void FaceSolver::ScanBcFaceDetail( IntSet& bcVertex, int bcType, int bcNameId )
{
    int nFaces = this->faceTopo->GetLeftCells().size();
    std::cout << " nFaces = " << nFaces << "\n";

    int nTraditionalBc = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int rCell = ( this->faceTopo->GetRightCells() )[ iFace ];

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
        int originalBcType = this->faceBcType[ iFace ];

        if ( originalBcType == ONEFLOW::INVALID_INDEX )
        {
            if ( this->CheckBcFace( bcVertex, ( this->faceTopo->GetFaces() )[ iFace ] ) )
            {
                ++ nBFaces;

                this->faceBcType[ iFace ] = bcType;
                this->faceBcKey[ iFace ] = bcNameId;
            }
        }
    }

    std::cout << " nFinalBcFace = " << nBFaces << " bcType = " << bcType << std::endl;

}

void FaceSolver::ScanInterfaceBc()
{
    int nFaces = this->faceTopo->GetLeftCells().size();

    int bcNameId = -1;
    int nInterFace = 0;
    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        int originalBcType = this->faceBcType[ iFace ];

        if ( originalBcType == ONEFLOW::INVALID_INDEX )
        {
            nInterFace ++;
            this->faceBcType[ iFace ] = BCTypeNull;
            this->faceBcKey[ iFace ] = bcNameId;
        }
    }

    std::cout << " nInterFace = " << nInterFace << std::endl;
}

int FaceSolver::GetNSimpleFace()
{
    int nSimpleFace = 0;

    for ( int iFace = 0; iFace < this->faceTopo->GetFaces().size(); ++ iFace )
    {
        int nCFace = this->childFid[ iFace ].size();
        if ( nCFace == 0 )
        {
            ++ nSimpleFace;
        }
    }
    return nSimpleFace;
}


EndNameSpace
