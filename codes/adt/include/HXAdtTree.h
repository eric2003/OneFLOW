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

#pragma once

#include "HXVector.h"
#include <memory>
#include <vector>
#include <cstring>
#include <algorithm>

#ifndef _WINDOWS
#include <string.h>
#endif

BeginNameSpace( ONEFLOW )

// ============================================================================
// Class Declaration & Implementation: HXAdtNode
// ============================================================================

template < typename T, typename U >
class HXAdtNode 
{
public:
    using AdtNode         = HXAdtNode<T, U>;
    using AdtNodeList     = HXVector<AdtNode*>;
    using AdtNodeListIter = typename AdtNodeList::iterator;

public:
    HXVector<U>  point;    // Coordinate array managed automatically by HXVector
    int          level;    // Depth level of the node in the tree
    AdtNode    * left;     // Pointer to the left child node (non-owning)
    AdtNode    * right;    // Pointer to the right child node (non-owning)
    T            item;     // User data storage (payload)
    int          dim;      // Spatial dimension count

public:
    explicit HXAdtNode( int dim = 3 )
        : level( 0 ), left( nullptr ), right( nullptr ), dim( dim )
    {
        point.resize( dim, static_cast<U>( 0 ) );
    }

    HXAdtNode( int dim, const U * coordinate, T data )
        : level( 0 ), left( nullptr ), right( nullptr ), item( std::move( data ) ), dim( dim )
    {
        point.resize( dim );
        std::memcpy( point.data(), coordinate, dim * sizeof( U ) );
    }

    // Centralized memory management via HXAdtTree avoids recursive destruction stack overflow
    ~HXAdtNode() = default;

    // Prevent accidental copying
    HXAdtNode( const HXAdtNode& ) = delete;
    HXAdtNode& operator=( const HXAdtNode& ) = delete;

    // Allow move semantics
    HXAdtNode( HXAdtNode&& ) noexcept = default;
    HXAdtNode& operator=( HXAdtNode&& ) noexcept = default;

    // Insert a child node into the current node
    void AddNode( AdtNode * node, U * nwmin, U * nwmax, const int & dim )
    {
        const int axis = level % dim;
        const U mid = static_cast<U>( 0.5 ) * ( nwmin[ axis ] + nwmax[ axis ] );

        if ( node->point[ axis ] <= mid )
        {
            if ( left != nullptr )
            {
                const U originalMax = nwmax[ axis ];
                nwmax[ axis ] = mid;
                left->AddNode( node, nwmin, nwmax, dim );
                nwmax[ axis ] = originalMax; // Backtrack
            }
            else
            {
                left = node;
                node->level = level + 1;
            }
        }
        else
        {
            if ( right != nullptr )
            {
                const U originalMin = nwmin[ axis ];
                nwmin[ axis ] = mid;
                right->AddNode( node, nwmin, nwmax, dim );
                nwmin[ axis ] = originalMin; // Backtrack
            }
            else
            {
                right = node;
                node->level = level + 1;
            }
        }
    }

    // Check if the current node is within the query bounding box
    [[nodiscard]] bool IsInRegion( const U * pmin, const U * pmax, const int & dim ) const
    {
        for ( int i = 0; i < dim; ++ i )
        {
            if ( point[ i ] < pmin[ i ] || point[ i ] > pmax[ i ] )
            {
                return false;
            }
        }
        return true;
    }

    // Retrieve all nodes located inside region (pmin, pmax)
    void FindNodesInRegion( const U * pmin, const U * pmax, U * nwmin, U * nwmax, const int & dim, AdtNodeList & ld ) const
    {
        if ( IsInRegion( pmin, pmax, dim ) )
        {
            ld.push_back( const_cast<AdtNode*>( this ) );
        }

        const int axis = level % dim;
        const U mid = static_cast<U>( 0.5 ) * ( nwmin[ axis ] + nwmax[ axis ] );

        if ( left != nullptr )
        {
            if ( pmin[ axis ] <= mid && pmax[ axis ] >= nwmin[ axis ] )
            {
                const U temp = nwmax[ axis ];
                nwmax[ axis ] = mid;
                left->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
                nwmax[ axis ] = temp; // Backtrack
            }
        }

        if ( right != nullptr )
        {
            if ( pmax[ axis ] >= mid && pmin[ axis ] <= nwmax[ axis ] )
            {
                const U temp = nwmin[ axis ];
                nwmin[ axis ] = mid;
                right->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
                nwmin[ axis ] = temp; // Backtrack
            }
        }
    }

    [[nodiscard]] int nCount() const
    {
        int iCount = 1;
        if ( this->left != nullptr )  iCount += left->nCount();
        if ( this->right != nullptr ) iCount += right->nCount();
        return iCount;
    }

    [[nodiscard]] T GetData() const { return item; }
};


// ============================================================================
// Class Declaration & Implementation: HXAdtTree
// ============================================================================

template < typename T, typename U >
class HXAdtTree
{
public:
    using AdtNode         = typename HXAdtNode<T, U>::AdtNode;
    using AdtNodeList     = typename HXAdtNode<T, U>::AdtNodeList;
    using AdtNodeListIter = typename HXAdtNode<T, U>::AdtNodeListIter;
    using AdtTree         = HXAdtTree<T, U>;

public:
    explicit HXAdtTree( int dim = 3 )
        : dim( dim ), root( nullptr )
    {
        pmin.resize( dim, static_cast<U>( 0.0 ) );
        pmax.resize( dim, static_cast<U>( 1.0 ) );
    }

    HXAdtTree( int dim, const U * pmin_in, const U * pmax_in )
        : dim( dim ), root( nullptr )
    {
        pmin.resize( dim );
        pmax.resize( dim );
        std::memcpy( pmin.data(), pmin_in, dim * sizeof( U ) );
        std::memcpy( pmax.data(), pmax_in, dim * sizeof( U ) );
    }

    HXAdtTree( int dim, const HXVector<U> & pmin_in, const HXVector<U> & pmax_in )
        : dim( dim ), pmin( pmin_in ), pmax( pmax_in ), root( nullptr )
    {
    }

    // Destructor automatically cleans up all nodes safely via unique_ptr
    ~HXAdtTree() = default;

    // Prevent accidental copying
    HXAdtTree( const HXAdtTree& ) = delete;
    HXAdtTree& operator=( const HXAdtTree& ) = delete;

    // Allow move semantics
    HXAdtTree( HXAdtTree&& ) noexcept = default;
    HXAdtTree& operator=( HXAdtTree&& ) noexcept = default;

    // Insert a node into the ADT tree (tree takes memory ownership)
    void AddNode( AdtNode * node )
    {
        if ( node == nullptr ) return;

        // Take ownership of memory to prevent memory leaks and stack overflows
        ownedNodes.emplace_back( node );

        if ( root == nullptr )
        {
            root = node;
            return;
        }

        // Allocate local vectors to ensure thread-safe backtracking without global side effects
        HXVector<U> localNwmin = this->pmin;
        HXVector<U> localNwmax = this->pmax;

        root->AddNode( node, localNwmin.data(), localNwmax.data(), dim );
    }

    // Find all nodes falling inside region (pmin, pmax)
    void FindNodesInRegion( const U * pmin_in, const U * pmax_in, AdtNodeList & ld ) const
    {
        if ( root == nullptr )
        {
            return;
        }

        HXVector<U> localNwmin = this->pmin;
        HXVector<U> localNwmax = this->pmax;

        root->FindNodesInRegion( pmin_in, pmax_in, localNwmin.data(), localNwmax.data(), dim, ld );
    }

    [[nodiscard]] int nCount() const
    {
        return root ? root->nCount() : 0;
    }

    // Get min bounding coordinates
    [[nodiscard]] U * GetMin() const
    {
        return const_cast<U*>( pmin.data() );
    }

    // Get max bounding coordinates
    [[nodiscard]] U * GetMax() const
    {
        return const_cast<U*>( pmax.data() );
    }

protected:
    int dim;
    HXVector<U> pmin; // Min boundary array
    HXVector<U> pmax; // Max boundary array
    AdtNode * root;

    // Centrally managed node memory: prevents memory leaks and recursive destruction stack overflows
    std::vector< std::unique_ptr<AdtNode> > ownedNodes; 
};

EndNameSpace