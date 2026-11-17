/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "InterfaceEndpoint.h"

#include "FaceTopo.h"
#include "InterFace.h"
#include "SolverDef.h"
#include "UnsGrid.h"

BeginNameSpace( ONEFLOW )

bool BuildInterfaceEndpoint(
    UnsGrid * grid,
    int neighborZone,
    int sendRecv,
    BoundaryEndpointRole role,
    BoundaryEndpointMap & endpoint )
{
    if ( grid == nullptr || grid->interFace == nullptr ) return false;

    InterFace * interFace = grid->interFace.get();
    int neighborIndex = -1;
    if ( ! interFace->TryGetNeighborIndex( neighborZone, neighborIndex ) )
        return false;

    const IntField & interfaceIds =
        interFace->GetInterfaceId( neighborIndex, sendRecv );
    endpoint.cellCount = grid->nCells + grid->nBFaces;
    endpoint.ownedCellCount = grid->nCells;
    endpoint.role = role;
    endpoint.cellIds.resize( interfaceIds.size() );

    for ( HXSize_t iFace = 0; iFace < interfaceIds.size(); ++ iFace )
    {
        int cellId = -1;
        if ( sendRecv == GREAT_SEND )
            grid->GetFaceTopo().GetSId( interfaceIds[ iFace ], 1, cellId );
        else
            grid->GetFaceTopo().GetTId( interfaceIds[ iFace ], 1, cellId );
        endpoint.cellIds[ iFace ] = cellId;
    }

    endpoint.Validate();
    return true;
}

EndNameSpace
