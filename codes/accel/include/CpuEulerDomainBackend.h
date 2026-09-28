#pragma once

#include "EulerDomain.h"

BeginNameSpace( ONEFLOW )

// Host-resident Euler domain state backend. Positive Advance is a host-side
// stage scheduler: the main solver remains authoritative for its MRField and
// supplies the numerical stage callback.
class CpuEulerDomainBackend final : public EulerDomainBackend
{
public:
    const char * Name() const override;
    bool IsAccelerator() const override;

    std::unique_ptr< EulerDomainState > CreateState(
        const EulerDomainProblem & problem,
        const EulerDomainStateKey & key ) const override;
    void Upload(
        EulerDomainState & state,
        const EulerDomainConstFieldView & field ) const override;
    void Advance(
        EulerDomainState & state,
        int steps,
        const EulerDomainRunOptions & options ) const override;
    void Download(
        const EulerDomainState & state,
        EulerDomainFieldView & field ) const override;
};

EndNameSpace
