#include "GridLayout.h"

#include <stdexcept>
#include <unordered_set>

BeginNameSpace( ONEFLOW )

void GridLayout::Validate() const
{
    std::unordered_set< int > pointIds;
    std::unordered_set< int > lineIds;
    std::unordered_set< int > faceIds;
    std::unordered_set< int > boundaryIds;

    for ( const auto & point : points )
    {
        if ( point.id <= 0 || ! pointIds.insert( point.id ).second )
        {
            throw std::runtime_error( "Invalid or duplicate grid point id" );
        }
    }

    for ( const auto & line : lines )
    {
        if ( line.id <= 0 || ! lineIds.insert( line.id ).second )
        {
            throw std::runtime_error( "Invalid or duplicate grid line id" );
        }
        if ( ! pointIds.contains( line.p1 ) || ! pointIds.contains( line.p2 ) )
        {
            throw std::runtime_error( "Grid line references an unknown point" );
        }
    }

    for ( const auto & circle : circles )
    {
        if ( circle.id <= 0 || ! lineIds.insert( circle.id ).second )
        {
            throw std::runtime_error( "Invalid or duplicate grid curve id" );
        }
        if ( ! pointIds.contains( circle.p1 ) ||
             ! pointIds.contains( circle.pc ) ||
             ! pointIds.contains( circle.p2 ) )
        {
            throw std::runtime_error( "Grid circle references an unknown point" );
        }
    }

    for ( const auto & dimension : dimensions )
    {
        if ( dimension.id <= 0 || ! lineIds.contains( dimension.id ) ||
             dimension.pointCount < 2 )
        {
            throw std::runtime_error( "Invalid grid dimension definition" );
        }
    }

    for ( const auto & distribution : distributions )
    {
        if ( ! lineIds.contains( distribution.lineId ) )
        {
            throw std::runtime_error( "Grid distribution references an unknown curve" );
        }
        if ( distribution.type == GridDistributionType::Copy )
        {
            for ( const int lineId : distribution.copyLineIds )
            {
                if ( ! lineIds.contains( lineId ) )
                {
                    throw std::runtime_error(
                        "Grid distribution copy references an unknown curve" );
                }
            }
        }
    }

    for ( const auto & boundary : boundaries )
    {
        if ( boundary.id <= 0 ||
             ! lineIds.contains( boundary.id ) ||
             ! boundaryIds.insert( boundary.id ).second )
        {
            throw std::runtime_error( "Invalid or duplicate grid boundary definition" );
        }
    }

    for ( const auto & relation : lineToFaces )
    {
        if ( relation.faceId <= 0 ||
             relation.position <= 0 ||
             ! lineIds.contains( relation.lineId ) )
        {
            throw std::runtime_error( "Invalid line-to-face relation" );
        }
        faceIds.insert( relation.faceId );
    }

    for ( const auto & relation : faceToBlocks )
    {
        if ( relation.blockId <= 0 ||
             relation.position <= 0 ||
             ! faceIds.contains( relation.faceId ) )
        {
            throw std::runtime_error( "Invalid face-to-block relation" );
        }
    }
}

EndNameSpace
