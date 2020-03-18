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


#pragma once
#include "HXVector.h"
using namespace std;

BeginNameSpace( ONEFLOW )

template < typename T, typename U >
class HXAdtNode 
{
public:
    typedef HXAdtNode< T, U >   AdtNode;
    typedef HXVector< AdtNode * > AdtNodeList;
    typedef typename AdtNodeList::iterator AdtNodeListIter;
public:
    U            * point;              // the coordinate of the node
    int          level;                // the level in the tree
    AdtNode      * left;               // the left tree
    AdtNode      * right;              // the right tree
    T            item;                 // any data stored
    int          dim;
public:
    HXAdtNode( int dim = 3 );
    HXAdtNode( int dim, U * coordinate, T data );
    ~HXAdtNode();

    // Add an Adt node under the current node 
    void AddNode( AdtNode * node, U * nwmin, U * nwmax, const int & dim );
    // is the current node inside region ( pmin, pmax )?
    bool IsInRegion( U * pmin, U * pmax, const int & dim );
    // ld carries all the nodes inside region ( pmin, pmax )
    void FindNodesInRegion( U * pmin, U * pmax, U * nwmin, U * nwmax, const int & dim, AdtNodeList & ld );
    int  nCount();
    T GetData() const { return item; };
};

template < typename T, typename U >
class HXAdtTree
{
public:
    typedef typename HXAdtNode<T, U>::AdtNode          AdtNode;
    typedef typename HXAdtNode<T, U>::AdtNodeList      AdtNodeList;
    typedef typename HXAdtNode<T, U>::AdtNodeListIter  AdtNodeListIter;
    typedef          HXAdtTree<T, U>                   AdtTree;
public:
    HXAdtTree( int dim = 3 );
    HXAdtTree( int dim, U * pmin, U * pmax );
    HXAdtTree( int dim, HXVector< U > & pmin, HXVector< U > &pmax );
    ~HXAdtTree();

    // Add an Adt node to the AdtTree 
    void AddNode( AdtNode * node );
    // Find All nodes inside the region ( pmin, pmax ) from the tree
    void FindNodesInRegion( U * pmin, U * pmax, AdtNodeList & ld );
    int  nCount()
    { 
        if ( root )
        {
            return root->nCount();
        }
        else
        {
            return 0;
        }
    };
    // Get the min coordinates of the tree
    U  * GetMin() const;
    // Get the max coordinates of the tree
    U  * GetMax() const;
protected:
    int     dim;
    U    * pmin, * pmax;
    AdtNode * root;
};

EndNameSpace

#include "AdtTree.hpp"