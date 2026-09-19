#pragma once

#include "EulerDomainStateRegistry.h"

BeginNameSpace( ONEFLOW )

// Centralizes create/upload ordering for solver-owned domain state.
// State enters the registry only after creation and host upload succeed.
class EulerDomainStateLifecycle
{
public:
    static EulerDomainState& Initialize(
        EulerDomainStateRegistry& registry,
        const EulerDomainBackend& backend,
        const EulerDomainProblem& problem,
        const EulerDomainStateKey& key,
        const EulerDomainConstFieldView& field );

    static EulerDomainState& Restart(
        EulerDomainStateRegistry& registry,
        const EulerDomainBackend& backend,
        const EulerDomainProblem& problem,
        const EulerDomainStateKey& key,
        const EulerDomainConstFieldView& field );
};

EndNameSpace
