#pragma once

#include "HXDefine.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

// Build per-physical-cell neighbor bounds as offsets from each cell value.
// The caller owns and clears the workspace to preserve the ghost-slot contract.
// The predicate keeps boundary classification outside this topology algorithm.
template< typename InterfaceBoundaryPredicate >
void BuildCellNeighborOffsetBounds(
    const RealField & q,
    const IntField & leftCells,
    const IntField & rightCells,
    const IntField & boundaryTypes,
    int nCells,
    int nBoundaryFaces,
    int nFaces,
    InterfaceBoundaryPredicate isInterfaceBoundary,
    RealField & minOffsets,
    RealField & maxOffsets,
    Real & minOffsetRange,
    Real & maxOffsetRange )
{
    for ( int cId = 0; cId < nCells; ++ cId )
    {
        minOffsets[ cId ] = q[ cId ];
        maxOffsets[ cId ] = q[ cId ];
    }

    // Interface boundary faces contribute the remote value to the local cell.
    for ( int fId = 0; fId < nBoundaryFaces; ++ fId )
    {
        if ( ! isInterfaceBoundary( boundaryTypes[ fId ] ) ) continue;

        const int lc = leftCells[ fId ];
        const int rc = rightCells[ fId ];
        minOffsets[ lc ] = MIN( minOffsets[ lc ], q[ rc ] );
        maxOffsets[ lc ] = MAX( maxOffsets[ lc ], q[ rc ] );
    }

    // Ordinary interior faces contribute to both adjacent physical cells.
    for ( int fId = nBoundaryFaces; fId < nFaces; ++ fId )
    {
        const int lc = leftCells[ fId ];
        const int rc = rightCells[ fId ];
        minOffsets[ lc ] = MIN( minOffsets[ lc ], q[ rc ] );
        maxOffsets[ lc ] = MAX( maxOffsets[ lc ], q[ rc ] );
        minOffsets[ rc ] = MIN( minOffsets[ rc ], q[ lc ] );
        maxOffsets[ rc ] = MAX( maxOffsets[ rc ], q[ lc ] );
    }

    // Only physical-cell entries are converted and included in the reduction.
    for ( int cId = 0; cId < nCells; ++ cId )
    {
        minOffsets[ cId ] -= q[ cId ];
        maxOffsets[ cId ] -= q[ cId ];
    }

    minOffsetRange = LARGE;
    maxOffsetRange = - LARGE;
    for ( int cId = 0; cId < nCells; ++ cId )
    {
        minOffsetRange = MIN( minOffsetRange, minOffsets[ cId ] );
        maxOffsetRange = MAX( maxOffsetRange, maxOffsets[ cId ] );
    }
}

EndNameSpace
