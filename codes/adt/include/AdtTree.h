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
#include <cstring>

#ifndef _WINDOWS
#include <string.h>
#endif

BeginNameSpace( ONEFLOW )

template < typename T, typename U >
class HXAdtNode 
{
public:
    using AdtNode       = HXAdtNode<T, U>;
    using AdtNodeList   = HXVector<AdtNode*>;
    using AdtNodeListIter = typename AdtNodeList::iterator;

public:
    HXVector<U>  point;            // Replaced raw pointer with HXVector for automatic memory management
    int          level;            // The level in the tree
    AdtNode    * left;             // The left child (non-owning raw pointer)
    AdtNode    * right;            // The right child (non-owning raw pointer)
    T            item;             // Any data stored
    int          dim;
public:
    HXAdtNode( int dim = 3 );
    HXAdtNode( int dim, U * coordinate, T data );

    // Intentionally does NOT delete left/right. 
    // Memory is centrally managed by HXAdtTree to prevent recursive stack overflow.
    ~HXAdtNode();

    // Add an Adt node under the current node 
    void AddNode( AdtNode * node, U * nwmin, U * nwmax, const int & dim );

    // Is the current node inside region ( pmin, pmax )?
    bool IsInRegion( U * pmin, U * pmax, const int & dim ) const;

    // ld carries all the nodes inside region ( pmin, pmax )
    void FindNodesInRegion( U * pmin, U * pmax, U * nwmin, U * nwmax, const int & dim, AdtNodeList & ld ) const;

    int nCount() const;
    T   GetData() const { return item; }
};


template < typename T, typename U >
class HXAdtTree
{
public:
    using AdtNode         = typename HXAdtNode<T, U>::AdtNode;
    using AdtNodeList     = typename HXAdtNode<T, U>::AdtNodeList;
    using AdtNodeListIter = typename HXAdtNode<T, U>::AdtNodeListIter;
    using AdtTree         = HXAdtTree<T, U>;

public:
    HXAdtTree( int dim = 3 );
    HXAdtTree( int dim, U * pmin_in, U * pmax_in );
    HXAdtTree( int dim, HXVector< U > & pmin_in, HXVector< U > & pmax_in );

    // Destructor relies on ownedNodes to safely clean up all nodes without recursion
    ~HXAdtTree();

    // Add an Adt node to the AdtTree (Tree takes ownership of the node)
    void AddNode( AdtNode * node );

    // Find all nodes inside the region ( pmin, pmax ) from the tree
    void FindNodesInRegion( U * pmin_in, U * pmax_in, AdtNodeList & ld ) const;

    int nCount() const;

    // Get the min coordinates of the tree
    U * GetMin() const;

    // Get the max coordinates of the tree
    U * GetMax() const;

protected:
    int dim;
    HXVector<U> pmin; // Replaced raw pointer with HXVector
    HXVector<U> pmax; // Replaced raw pointer with HXVector
    AdtNode * root;

    // Centralized memory management: owns all nodes, preventing memory leaks and recursive destruction issues
    std::vector< std::unique_ptr<AdtNode> > ownedNodes; 
};


EndNameSpace

#include "AdtTree.hpp"
