/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include "BoundaryExchangeContract.h"

BeginNameSpace( ONEFLOW )

class UnsGrid;

// Resolve one zone interface into the cell endpoint used by host-staged or
// device-direct exchange.  The topology mapping is shared by all backends;
// only the transfer implementation is backend-specific.
bool BuildInterfaceEndpoint(
    UnsGrid * grid,
    int neighborZone,
    int sendRecv,
    BoundaryEndpointRole role,
    BoundaryEndpointMap & endpoint );

EndNameSpace
