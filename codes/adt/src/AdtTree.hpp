/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2025 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

// ============================================================================
// Template Implementations: HXAdtNode
// ============================================================================

template < typename T, typename U >
HXAdtNode<T,U>::HXAdtNode( int dim )
{
    this->dim = dim;
    point.resize(dim, static_cast<U>(0));
    level = 0;
    left  = nullptr;
    right = nullptr;
}

template < typename T, typename U >
HXAdtNode<T,U>::HXAdtNode( int dim, U * coordinate, T data )
{
    this->dim = dim;
    point.resize(dim);
    std::memcpy( point.data(), coordinate, dim * sizeof( U ) );

    level = 0;
    left  = nullptr;
    right = nullptr;
    item  = data;
}

template < typename T, typename U >
HXAdtNode<T,U>::~HXAdtNode()
{
    // Do NOT delete left or right here. 
    // HXAdtTree centrally manages memory via std::unique_ptr to avoid deep recursion stack overflow.
    // HXVector 'point' cleans up its own memory automatically.
}

template < typename T, typename U >
int HXAdtNode<T,U>::nCount() const
{
    int iCount = 1;
    if ( this->left != nullptr )
    {
        iCount += left->nCount();
    }
    if ( this->right != nullptr )
    {
        iCount += right->nCount();
    }
    return iCount;
}

// Add an Adt node under the current node
template < typename T, typename U >
void HXAdtNode<T,U>::AddNode( AdtNode * node, U * nwmin, U * nwmax, const int & dim )
{
    int axis = level % dim;
    U mid = static_cast<U>(0.5) * ( nwmin[ axis ] + nwmax[ axis ] );

    if ( node->point[ axis ] <= mid )
    {
        if ( left != nullptr )
        {
            U originalMax = nwmax[ axis ];
            nwmax[ axis ] = mid;
            left->AddNode( node, nwmin, nwmax, dim );
            nwmax[ axis ] = originalMax; // Backtrack
        }
        else
        {
            left        = node;
            node->level = level + 1;
        }
    }
    else
    {
        if ( right != nullptr )
        {
            U originalMin = nwmin[ axis ];
            nwmin[ axis ] = mid;
            right->AddNode( node, nwmin, nwmax, dim );
            nwmin[ axis ] = originalMin; // Backtrack
        }
        else
        {
            right       = node;
            node->level = level + 1;
        }
    }
}

// Is the current node inside region ( pmin, pmax )?
template < typename T, typename U >
bool HXAdtNode<T,U>::IsInRegion( U * pmin, U * pmax, const int & dim ) const
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
void HXAdtNode<T,U>::FindNodesInRegion( U * pmin, U * pmax, U * nwmin, U * nwmax, const int & dim, AdtNodeList & ld ) const
{
    if ( IsInRegion( pmin, pmax, dim ) )
    {
        // Cast to non-const pointer to match the original AdtNodeList type requirement
        ld.push_back( const_cast<AdtNode*>(this) );
    }

    int axis = level % dim;
    U mid = static_cast<U>(0.5) * ( nwmin[ axis ] + nwmax[ axis ] );

    if ( left != nullptr )
    {
        if ( pmin[ axis ] <= mid && pmax[ axis ] >= nwmin[ axis ] )
        {
            U temp = nwmax[ axis ];
            nwmax[ axis ] = mid;
            left->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
            nwmax[ axis ] = temp; // Backtrack
        }
    }

    if ( right != nullptr )
    {
        if ( pmax[ axis ] >= mid && pmin[ axis ] <= nwmax[ axis ] )
        {
            U temp = nwmin[ axis ];
            nwmin[ axis ] = mid;
            right->FindNodesInRegion( pmin, pmax, nwmin, nwmax, dim, ld );
            nwmin[ axis ] = temp; // Backtrack
        }
    }
}

// ============================================================================
// Template Implementations: HXAdtTree
// ============================================================================

template < typename T, typename U >
HXAdtTree<T,U>::HXAdtTree( int dim )
{
    this->dim = dim;
    pmin.resize(dim, static_cast<U>(0.0));
    pmax.resize(dim, static_cast<U>(1.0));
    root = nullptr;
}

template < typename T, typename U >
HXAdtTree<T,U>::HXAdtTree( int dim, U * pmin_in, U * pmax_in )
{
    this->dim = dim;
    pmin.resize(dim);
    pmax.resize(dim);
    for ( int i = 0; i < dim; ++ i )
    { 
        this->pmin[ i ] = pmin_in[ i ];
        this->pmax[ i ] = pmax_in[ i ];
    }
    root = nullptr;
}

template < typename T, typename U >
HXAdtTree<T,U>::HXAdtTree( int dim, HXVector< U > & pmin_in, HXVector< U > & pmax_in )
{
    this->dim = dim;
    pmin = pmin_in;
    pmax = pmax_in;
    root = nullptr;
}

template < typename T, typename U >
HXAdtTree<T,U>::~HXAdtTree()
{  
    // std::vector of unique_ptr automatically cleans up all owned nodes safely.
    // No manual delete[] or recursive delete is needed, preventing stack overflow.
}

// Add an Adt node to the AdtTree 
template < typename T, typename U >
void HXAdtTree<T,U>::AddNode( AdtNode * node )
{
    // Take ownership of the raw pointer immediately to prevent memory leaks
    ownedNodes.emplace_back( node );

    if ( root == nullptr )
    {
        root = node;
        return;
    }

    // Use local HXVector to avoid heap allocation (new/delete) during recursion
    HXVector<U> localNwmin = this->pmin;
    HXVector<U> localNwmax = this->pmax;

    root->AddNode( node, localNwmin.data(), localNwmax.data(), dim );
}

// Find all nodes inside the region ( pmin, pmax ) from the tree
template < typename T, typename U >
void HXAdtTree<T,U>::FindNodesInRegion( U * pmin_in, U * pmax_in, AdtNodeList & ld ) const
{
    if ( root == nullptr )
    {
        return;
    }

    // Use local HXVector to avoid heap allocation (new/delete) during recursion
    HXVector<U> localNwmin = this->pmin;
    HXVector<U> localNwmax = this->pmax;

    root->FindNodesInRegion( pmin_in, pmax_in, localNwmin.data(), localNwmax.data(), dim, ld );
}

template < typename T, typename U >
int HXAdtTree<T,U>::nCount() const
{ 
    if ( root != nullptr )
    {
        return root->nCount();
    }
    return 0;
}

template < typename T, typename U >
U * HXAdtTree<T,U>::GetMin() const
{
    // Return pointer to internal data. 
    // Note: Caller should not modify this data, but signature is kept for backward compatibility.
    return const_cast<U*>(pmin.data());
}

template < typename T, typename U >
U * HXAdtTree<T,U>::GetMax() const
{
    // Return pointer to internal data.
    return const_cast<U*>(pmax.data());
}


EndNameSpace
