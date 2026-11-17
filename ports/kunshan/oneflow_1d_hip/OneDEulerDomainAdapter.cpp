#include "OneDEulerDomainAdapter.h"

#include <cstring>
#include <stdexcept>
#include <utility>

namespace oneflow_1d
{
namespace
{

EulerBoundary MapBoundary( ONEFLOW::EulerDomainBoundary boundary )
{
    switch ( boundary )
    {
    case ONEFLOW::EulerDomainBoundary::Periodic:
        return EulerBoundary::Periodic;
    case ONEFLOW::EulerDomainBoundary::Outflow:
        return EulerBoundary::Transmissive;
    case ONEFLOW::EulerDomainBoundary::Wall:
    case ONEFLOW::EulerDomainBoundary::Inflow:
        break;
    }
    throw std::invalid_argument(
        "1D Euler adapter does not support the requested boundary contract" );
}

EulerProblem MapProblem( const ONEFLOW::EulerDomainProblem & problem )
{
    ONEFLOW::ValidateEulerDomainProblem( problem );
    if ( problem.nEquations != EulerComponents || problem.nGhostCells != 0 )
    {
        throw std::invalid_argument(
            "1D Euler adapter requires three equations and no external ghosts" );
    }
    return EulerProblem{
        problem.nCells,
        problem.gamma,
        problem.dt,
        problem.dx,
        MapBoundary( problem.boundary ),
        EulerMethod::Rusanov };
}

void ValidateKey(
    const ONEFLOW::EulerDomainStateKey & key,
    ONEFLOW::AccelBackendKind backendKind,
    int deviceId )
{
    if ( key.solverIndex < 0 || key.localZoneId < 0 || key.gridLevel < 0
         || key.backend != backendKind || key.deviceId != deviceId )
    {
        throw std::invalid_argument(
            "invalid 1D Euler adapter state key" );
    }
}

struct AdapterState final : ONEFLOW::EulerDomainState
{
    const EulerBackend * owner = nullptr;
    ONEFLOW::EulerDomainProblem problem;
    ONEFLOW::EulerDomainStateKey key;
    std::unique_ptr< EulerState > state;
};

AdapterState & AsAdapterState(
    ONEFLOW::EulerDomainState & state,
    const EulerBackend * owner,
    ONEFLOW::AccelBackendKind backendKind,
    int deviceId )
{
    auto * result = dynamic_cast< AdapterState * >( &state );
    if ( result == nullptr || result->owner != owner
         || result->key.backend != backendKind
         || result->key.deviceId != deviceId )
    {
        throw std::invalid_argument(
            "Euler state does not belong to this 1D domain adapter" );
    }
    return *result;
}

const AdapterState & AsAdapterState(
    const ONEFLOW::EulerDomainState & state,
    const EulerBackend * owner,
    ONEFLOW::AccelBackendKind backendKind,
    int deviceId )
{
    auto * result = dynamic_cast< const AdapterState * >( &state );
    if ( result == nullptr || result->owner != owner
         || result->key.backend != backendKind
         || result->key.deviceId != deviceId )
    {
        throw std::invalid_argument(
            "Euler state does not belong to this 1D domain adapter" );
    }
    return *result;
}

EulerRunOptions MapRunOptions(
    int steps,
    const ONEFLOW::EulerDomainRunOptions & options )
{
    if ( steps <= 0 )
    {
        throw std::invalid_argument(
            "1D Euler adapter advance steps must be positive" );
    }
    if ( options.stageCount != 1 || options.stageCallback != nullptr )
    {
        throw std::invalid_argument(
            "1D Euler adapter owns its internal RK stages" );
    }

    EulerRunOptions result;
    if ( options.mode == ONEFLOW::EulerDomainRunMode::FullTrace )
    {
        result.mode = EulerRunMode::FullTrace;
        result.trace = static_cast< EulerTrace * >( options.trace );
    }
    else
    {
        result.mode = EulerRunMode::NoTrace;
        if ( options.trace != nullptr )
        {
            throw std::invalid_argument(
                "NoTrace cannot provide a 1D Euler trace" );
        }
    }
    result.stats = static_cast< EulerRunStats * >( options.stats );
    return result;
}

}

EulerDomainBackendAdapter::EulerDomainBackendAdapter(
    const EulerBackend & backend,
    ONEFLOW::AccelBackendKind backendKind,
    int deviceId )
    : backend_( &backend ), backendKind_( backendKind ), deviceId_( deviceId )
{
    const bool isCpu = std::strcmp( backend.Name(), "CPU" ) == 0;
    const bool isHip = std::strcmp( backend.Name(), "HIP" ) == 0;
    if ( ( backendKind == ONEFLOW::AccelBackendKind::CPU && ! isCpu )
         || ( backendKind == ONEFLOW::AccelBackendKind::HIP && ! isHip )
         || ( backendKind != ONEFLOW::AccelBackendKind::CPU
              && backendKind != ONEFLOW::AccelBackendKind::HIP ) )
    {
        throw std::invalid_argument(
            "1D Euler adapter does not support the requested backend kind" );
    }
    if ( backend.IsAccelerator()
         != ( backendKind != ONEFLOW::AccelBackendKind::CPU ) )
    {
        throw std::invalid_argument(
            "1D Euler adapter backend kind does not match the wrapped backend" );
    }
    if ( ( backendKind == ONEFLOW::AccelBackendKind::CPU && deviceId != -1 )
         || ( backendKind != ONEFLOW::AccelBackendKind::CPU && deviceId < 0 ) )
    {
        throw std::invalid_argument(
            "1D Euler adapter has an invalid device identity" );
    }
}

const char * EulerDomainBackendAdapter::Name() const
{
    return backend_->Name();
}

bool EulerDomainBackendAdapter::IsAccelerator() const
{
    return backend_->IsAccelerator();
}

std::unique_ptr< ONEFLOW::EulerDomainState >
EulerDomainBackendAdapter::CreateState(
    const ONEFLOW::EulerDomainProblem & problem,
    const ONEFLOW::EulerDomainStateKey & key ) const
{
    const EulerProblem portProblem = MapProblem( problem );
    ValidateKey( key, backendKind_, deviceId_ );

    auto result = std::make_unique< AdapterState >();
    result->owner = backend_;
    result->problem = problem;
    result->key = key;
    result->state = backend_->CreateState( portProblem );
    return result;
}

void EulerDomainBackendAdapter::Upload(
    ONEFLOW::EulerDomainState & state,
    const ONEFLOW::EulerDomainConstFieldView & field ) const
{
    AdapterState & adapter =
        AsAdapterState( state, backend_, backendKind_, deviceId_ );
    ONEFLOW::ValidateEulerDomainField( adapter.problem, field );
    backend_->Upload( *adapter.state, field.values );
}

void EulerDomainBackendAdapter::Advance(
    ONEFLOW::EulerDomainState & state,
    int steps,
    const ONEFLOW::EulerDomainRunOptions & options ) const
{
    AdapterState & adapter =
        AsAdapterState( state, backend_, backendKind_, deviceId_ );
    backend_->Advance(
        *adapter.state, steps, MapRunOptions( steps, options ) );
}

void EulerDomainBackendAdapter::Download(
    const ONEFLOW::EulerDomainState & state,
    ONEFLOW::EulerDomainFieldView & field ) const
{
    const AdapterState & adapter =
        AsAdapterState( state, backend_, backendKind_, deviceId_ );
    ONEFLOW::ValidateEulerDomainField( adapter.problem, field );
    backend_->Download( *adapter.state, field.values );
}

} // namespace oneflow_1d
