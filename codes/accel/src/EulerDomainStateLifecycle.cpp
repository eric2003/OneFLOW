#include "EulerDomainStateLifecycle.h"

#include <stdexcept>

BeginNameSpace( ONEFLOW )

namespace
{

EulerDomainState& CreateUploadAndRegister(
    EulerDomainStateRegistry& registry,
    const EulerDomainBackend& backend,
    const EulerDomainProblem& problem,
    const EulerDomainStateKey& key,
    const EulerDomainConstFieldView& field )
{
    ValidateEulerDomainField( problem, field );

    std::unique_ptr< EulerDomainState > state =
        backend.CreateState( problem, key );
    if ( state == nullptr )
    {
        throw std::runtime_error(
            "Euler domain backend returned a null state" );
    }

    backend.Upload( *state, field );
    registry.Insert( key, std::move( state ) );
    return registry.Get( key );
}

}

EulerDomainState& EulerDomainStateLifecycle::Initialize(
    EulerDomainStateRegistry& registry,
    const EulerDomainBackend& backend,
    const EulerDomainProblem& problem,
    const EulerDomainStateKey& key,
    const EulerDomainConstFieldView& field )
{
    registry.Invalidate( key );
    return CreateUploadAndRegister( registry, backend, problem, key, field );
}

EulerDomainState& EulerDomainStateLifecycle::Restart(
    EulerDomainStateRegistry& registry,
    const EulerDomainBackend& backend,
    const EulerDomainProblem& problem,
    const EulerDomainStateKey& key,
    const EulerDomainConstFieldView& field )
{
    // Keep restart explicit so a future device backend can observe the
    // invalidation boundary independently from initialization.
    registry.Invalidate( key );
    return CreateUploadAndRegister( registry, backend, problem, key, field );
}

EndNameSpace
