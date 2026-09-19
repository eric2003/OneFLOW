#pragma once

#include "EulerDomain.h"

BeginNameSpace( ONEFLOW )

class CpuEulerDomainBackend;
class SimuContext;

void SyncCurrentEulerDomainState(
    SimuContext & context,
    CpuEulerDomainBackend & backend,
    bool restart );

void UploadCurrentEulerDomainState(
    SimuContext & context,
    CpuEulerDomainBackend & backend );

void SyncAllEulerDomainStates( SimuContext & context );

EndNameSpace
