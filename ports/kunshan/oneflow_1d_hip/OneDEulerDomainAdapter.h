#pragma once

#include "EulerDomain.h"
#include "OneDEulerBackend.h"

namespace oneflow_1d
{

// Thin compatibility layer that exposes the validated 1D Euler backend
// through the shared solver lifecycle contract. The wrapped backend remains
// the numerical oracle; this adapter owns no second copy of CPU or device
// state.
class EulerDomainBackendAdapter final : public ONEFLOW::EulerDomainBackend
{
public:
    EulerDomainBackendAdapter(
        const EulerBackend & backend,
        ONEFLOW::AccelBackendKind backendKind,
        int deviceId );

    const char * Name() const override;
    bool IsAccelerator() const override;

    std::unique_ptr< ONEFLOW::EulerDomainState > CreateState(
        const ONEFLOW::EulerDomainProblem & problem,
        const ONEFLOW::EulerDomainStateKey & key ) const override;
    void Upload(
        ONEFLOW::EulerDomainState & state,
        const ONEFLOW::EulerDomainConstFieldView & field ) const override;
    void Advance(
        ONEFLOW::EulerDomainState & state,
        int steps,
        const ONEFLOW::EulerDomainRunOptions & options ) const override;
    void Download(
        const ONEFLOW::EulerDomainState & state,
        ONEFLOW::EulerDomainFieldView & field ) const override;

private:
    const EulerBackend * backend_ = nullptr;
    ONEFLOW::AccelBackendKind backendKind_ = ONEFLOW::AccelBackendKind::CPU;
    int deviceId_ = -1;
};

} // namespace oneflow_1d
