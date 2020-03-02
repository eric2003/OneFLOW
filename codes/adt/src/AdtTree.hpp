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
#include <memory>
#ifndef _WINDOWS
   #include <string.h>
#endif

using namespace std;

BeginNameSpace( ONEFLOW )

template < typename T, typename U >
HXAdtNode<T,U>::HXAdtNode( int dim )
{
    this->dim = dim;
    point = new U[ dim ];
    level = 0;
    left  = 0;
    right = 0;
}
  
template < typename T, typename U >
HXAdtNode<T,U>::HXAdtNode( int dim, U * coordinate, T data )
{
    this->dim = dim;
    point = new U [ dim ];
    memcpy( point, coordinate, dim * sizeof( U ) );

    level = 0   ;
    left  = 0   ;
    right = 0   ;
    item  = data;
}

template < typename T, typename U >
HXAdtNode<T,U>::~HXAdtNode()
{
    delete [] point;
    delete left;
    delete right;
}

template < typename T, typename U >
int HXAdtNode<T,U>::nCount()
{
    int iCount = 0;
    iCount += 1;
    if ( this->left )
    {
        iCount += left->nCount();
    }
    if ( this->right )
    {
        iCount += right->nCount();
    }
    return iCount;
}

// Add an Adt node under the current node
template < typename T, typename U >
void HXAdtNode<T,U>::AddNode( AdtNode * node, U * nwmin, U * nwmax, const int & dim )
{
    int axis  = level % dim;
    U mid = 0.5 * ( nwmin[ axis ] + nwmax[ axis ] );
      
    if ( node->point[ axis ] <= mid )
    {
        if ( left )
        {
            nwmax[ axis ] = mid;
            left->AddNode( node, nwmin, nwmax, dim );
        }
        else
        {
            left        = node     ;
            node->level = level + 1;
        }
    }
    else
    {
        if ( right )
        {
            nwmin[ axis ] = mid;
            right->AddNode( node, nwmin, nwmax, dim );
        }
        else
        {
            right       = node     ;
            node->level = level + 1;
        }
    }
}

// is the current node inside region ( pmin, pmax )?
template < typename T, typename U >
bool HXAdtNode<T,U>::IsInRegion( U * pmin, U * pmax, const int & dim )
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

// ld carries all the nodes inside region ( pmin, pmax )
template < typename T, typename U >
void HXAdtNode<T,U>::FindNodesInRegion( U * pmin, U * pmax, U * nwmin, U * nwmax, const int & dim, AdtNodeList & ld )
{
    int     axis;
    U       mid, temp;

    if ( IsInRegion( pmin, pmax, dim ) )
    {
        ld.push_back( this );
    }

    axis = level%dim;
    mid = 0.5 * ( nwmin[ axis ] + nwmax[ axis ] );
      
    if ( left )
    {
        if ( pmin[ axis ] <= mid && pmax[ axis ] >= nwmin[ axis ] )
        {
            temp        = nwmax[ axis ];
            nwmax[ axis ] = mid;
            left->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
            nwmax[ axis ] = temp;
        }
    }

    if ( right )
    {
        if ( pmax[ axis ] >= mid && pmin[ axis ] <= nwmax[ axis ] )
        {
            temp        = nwmin[ axis ];
            nwmin[ axis ] = mid;
            right->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
            nwmin[ axis ] = temp;
        }
    }

    return;
}


template < typename T, typename U >
HXAdtTree<T,U>::HXAdtTree( int dim )
{
    this->dim = dim;
    pmin = new U[ dim ];
    pmax = new U[ dim ];
    for ( int i = 0; i < dim; ++ i )
    { 
        pmin[ i ] = 0.0;
        pmax[ i ] = 1.0;
    }
    root = 0;
}

template < typename T, typename U >
HXAdtTree<T,U>::HXAdtTree( int dim, U * pmin, U * pmax )
{
    this->dim = dim;
    this->pmin = new U[ dim ];
    this->pmax = new U[ dim ];
    for ( int i = 0; i < dim; ++ i )
    { 
        this->pmin[ i ] = pmin[ i ];
        this->pmax[ i ] = pmax[ i ];
    }
    root = 0;
}

template < typename T, typename U >
HXAdtTree<T,U>::HXAdtTree( int dim, HXVector< U > & pmin, HXVector< U > & pmax )
{
    this->dim = dim;
    this->pmin = new U[ dim ];
    this->pmax = new U[ dim ];
    for ( int i = 0; i < dim; ++ i )
    { 
        this->pmin[ i ] = pmin[ i ];
        this->pmax[ i ] = pmax[ i ];
    }
    root = 0;
}

template < typename T, typename U >
HXAdtTree<T,U>::~HXAdtTree()
{  
    delete [] pmin;
    delete [] pmax;
    delete root;
}

// Add an Adt node to the AdtTree 
template < typename T, typename U >
void HXAdtTree<T,U>::AddNode( AdtNode * node )
{
    U * nwmin = new U [ dim ];
    U * nwmax = new U [ dim ];
    memcpy( nwmin, this->pmin, dim * sizeof( U ) );
    memcpy( nwmax, this->pmax, dim * sizeof( U ) );
      
    if ( root == 0 )
    {
        root = node;
    }
    else
    {
        root->AddNode( node, nwmin, nwmax, dim );
    }
        
    delete [] nwmin;
    delete [] nwmax;
}

// Find All nodes inside the region ( pmin, pmax ) from the tree
template < typename T, typename U >
void HXAdtTree<T,U>::FindNodesInRegion( U * pmin, U * pmax, AdtNodeList & ld )
{
    U * nwmin = new U [ dim ];
    U * nwmax = new U [ dim ];
    memcpy( nwmin, this->pmin, dim * sizeof( U ) );
    memcpy( nwmax, this->pmax, dim * sizeof( U ) );

    if ( root )
    {
        root->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
    }
    delete [] nwmin;
    delete [] nwmax;
}


// Get the min coordinates of the tree
template < typename T, typename U >
U * HXAdtTree<T,U>::GetMin() const
{
    return pmin;
}

// Get the max coordinates of the tree
template < typename T, typename U >
U * HXAdtTree<T,U>::GetMax() const
{
    return pmax;
}

EndNameSpace