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
#include "HXDefine.h"
#include "HXArray.h"
#include <memory>

BeginNameSpace( ONEFLOW )

// Cell-centered <-> node-centered conversion (logistics for multi-method coupling).
// Ownership: returned unique_ptr is owned by the caller; cell sources are non-owning reads.

std::unique_ptr<MRField> AllocateNodeField( int nEqu = 1 );

// Look up cell field by name (const read) and allocate a node field.
std::unique_ptr<MRField> InterpolateCellToNode( const std::string & name );

// Interpolate from an existing cell-centered component (const read).
std::unique_ptr<MRField> InterpolateCellToNode( const RealField & qc );

// Write node field from const cell component.
void InterpolateCellToNodeForComponent( RealField & qNodeField, const RealField & qField );

// BC correction on node accumulation; qField is read-only cell data.
void ApplyBoundaryConditionToNodeField( RealField & qNodeField, const RealField & qField, RealField & nCount, int bcType, bool twoSide );

// Symmetric logistics: node-centered -> cell-centered (caller owns result).
std::unique_ptr<MRField> AllocateCellField( int nEqu = 1 );
std::unique_ptr<MRField> InterpolateNodeToCell( const RealField & qn );
void InterpolateNodeToCellForComponent( RealField & qCellField, const RealField & qNodeField );

template < typename T >
void ReorderList( HXVector< T > & x, IntField & indexList )
{
    HXVector< T > tmp = x;
    for ( HXSize_t i = 0; i < indexList.size(); ++ i )
    {
        x[ i ] = tmp[ indexList[ i ] ];
    }
}

EndNameSpace
