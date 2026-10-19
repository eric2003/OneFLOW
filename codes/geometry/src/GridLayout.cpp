#include "GridLayout.h"

#include <set>
#include <stdexcept>
#include <utility>
#include <unordered_set>

BeginNameSpace( ONEFLOW )

GridBoundaryKind GridLayout::ClassifyBoundaryType( int boundaryType )
{
    if ( boundaryType == 0 ) return GridBoundaryKind::None;
    if ( boundaryType == -1 ) return GridBoundaryKind::Interface;
    if ( boundaryType == 100 ) return GridBoundaryKind::Interior;
    if ( boundaryType == 1000 ) return GridBoundaryKind::Overset;
    return GridBoundaryKind::Physical;
}

void GridLayout::Validate() const
{
    std::unordered_set< int > pointIds;
    std::unordered_set< int > lineIds;
    std::unordered_set< int > faceIds;
    std::unordered_set< int > boundaryIds;
    std::unordered_set< int > dimensionIds;
    std::set< std::pair< int, int > > lineToFaceSlots;
    std::set< std::pair< int, int > > blockFaceSlots;

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
             dimension.pointCount < 2 ||
             ! dimensionIds.insert( dimension.id ).second )
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
        else if ( distribution.type == GridDistributionType::BoundaryLayer )
        {
            if ( distribution.nearWallPointCount < 2 ||
                 distribution.nearWallDistance <= 0.0 )
            {
                throw std::runtime_error(
                    "Invalid boundary-layer distribution parameters" );
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
        if ( ! lineToFaceSlots.insert(
                  { relation.faceId, relation.position } ).second )
        {
            throw std::runtime_error( "Duplicate line-to-face position" );
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
        if ( ! blockFaceSlots.insert(
                  { relation.blockId, relation.position } ).second )
        {
            throw std::runtime_error( "Duplicate face-to-block position" );
        }
    }
}

EndNameSpace
