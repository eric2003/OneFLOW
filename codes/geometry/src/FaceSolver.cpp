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
#include <limits>
#include <stdexcept>
#include <unordered_set>


BeginNameSpace( ONEFLOW )

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

std::unique_ptr< FaceTopo > FaceSolver::ReleaseFaceTopo() noexcept
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

void FaceSolver::ScanPolygonFace( CgnsSection & cgnsSection )
{
    const auto & offsets = cgnsSection.ePosList;
    const auto & connectivity = cgnsSection.connList;

    if ( cgnsSection.nElement < 0 ||
         static_cast< std::size_t >( cgnsSection.nElement ) + 1 > offsets.size() )
    {
        throw std::runtime_error( "FaceSolver::ScanPolygonFace: polygon offsets are incomplete" );
    }

    // Validate every span before mutating the face lookup or topology.
    for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        if ( start < 0 || end < start ||
             static_cast< std::size_t >( end ) > connectivity.size() )
        {
            throw std::runtime_error( "FaceSolver::ScanPolygonFace: polygon connectivity span is invalid" );
        }
    }

    IntField faceNodes;
    for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        faceNodes.resize( 0 );
        for ( CgInt i = start; i < end; ++ i )
        {
            faceNodes.push_back( connectivity[ static_cast< std::size_t >( i ) ] );
        }

        auto [faceIndex, isNew] = faceLookup.FindOrAdd( faceNodes );
        if ( isNew )
        {
            this->faceTopo->GetFaces().push_back( faceNodes );
            this->faceTopo->GetFaceTypes().push_back( cgnsSection.eType );
            this->faceTopo->GetFaceFlags().push_back( 0 );
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

void FaceSolver::ScanPolyhedronElement( CgnsSection & cgnsSection )
{
    const auto & offsets = cgnsSection.ePosList;
    const auto & connectivity = cgnsSection.connList;
    const auto & faces = this->faceTopo->GetFaces();

    if ( cgnsSection.nElement < 0 ||
         static_cast< std::size_t >( cgnsSection.nElement ) + 1 > offsets.size() )
    {
        throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: element offsets are incomplete" );
    }

    // Validate the complete section before changing cell adjacency.
    for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        if ( start < 0 || end < start ||
             static_cast< std::size_t >( end ) > connectivity.size() )
        {
            throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: connectivity span is invalid" );
        }

        std::unordered_set< std::size_t > referencedFaces;
        for ( CgInt i = start; i < end; ++ i )
        {
            const CgInt signedFaceId = connectivity[ static_cast< std::size_t >( i ) ];
            if ( signedFaceId == std::numeric_limits< CgInt >::min() )
            {
                throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: face reference is out of range" );
            }

            const CgInt faceId = signedFaceId < 0 ? -signedFaceId : signedFaceId;
            if ( faceId < 0 || static_cast< std::size_t >( faceId ) >= faces.size() )
            {
                throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: face reference is out of range" );
            }

            if ( ! referencedFaces.insert( static_cast< std::size_t >( faceId ) ).second )
            {
                throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: polyhedron references the same face more than once" );
            }
        }
    }

    for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        if ( start < 0 || end < start ||
             static_cast< std::size_t >( end ) > connectivity.size() )
        {
            throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: connectivity span is invalid" );
        }

        for ( CgInt i = start; i < end; ++ i )
        {
            const CgInt signedFaceId = connectivity[ static_cast< std::size_t >( i ) ];
            if ( signedFaceId == std::numeric_limits< CgInt >::min() )
            {
                throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: face reference is out of range" );
            }

            const CgInt faceId = signedFaceId < 0 ? -signedFaceId : signedFaceId;
            if ( faceId < 0 || static_cast< std::size_t >( faceId ) >= faces.size() )
            {
                throw std::runtime_error( "FaceSolver::ScanPolyhedronElement: face reference is out of range" );
            }

            const std::size_t polygonFaceId = static_cast< std::size_t >( faceId );
            const int faceFlags = this->faceTopo->GetFaceFlags()[ polygonFaceId ];

            if ( faceFlags == 0 )
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
