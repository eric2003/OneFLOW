#include "CpuEulerDomainBackend.h"

#include <algorithm>
#include <stdexcept>
#include <utility>
#include <vector>

BeginNameSpace( ONEFLOW )

namespace
{

struct CpuEulerDomainState final : EulerDomainState
{
    EulerDomainProblem problem;
    EulerDomainStateKey key;
    std::vector< Real > values;
};

void ValidateKey( const EulerDomainStateKey & key )
{
    if ( key.solverIndex < 0 || key.localZoneId < 0 || key.gridLevel < 0
         || key.backend != AccelBackendKind::CPU )
    {
        throw std::invalid_argument(
            "invalid CPU Euler domain state key" );
    }
}

CpuEulerDomainState & AsCpuState( EulerDomainState & state )
{
    auto * cpuState = dynamic_cast< CpuEulerDomainState * >( &state );
    if ( cpuState == nullptr )
    {
        throw std::invalid_argument(
            "Euler state does not belong to the CPU domain backend" );
    }
    return *cpuState;
}

const CpuEulerDomainState & AsCpuState( const EulerDomainState & state )
{
    auto * cpuState = dynamic_cast< const CpuEulerDomainState * >( &state );
    if ( cpuState == nullptr )
    {
        throw std::invalid_argument(
            "Euler state does not belong to the CPU domain backend" );
    }
    return *cpuState;
}

std::size_t ValueCount( const EulerDomainProblem & problem )
{
    return static_cast< std::size_t >( problem.nCells )
        * static_cast< std::size_t >( problem.nEquations );
}

}

const char * CpuEulerDomainBackend::Name() const
{
    return "CPU-EulerDomain";
}

bool CpuEulerDomainBackend::IsAccelerator() const
{
    return false;
}

std::unique_ptr< EulerDomainState > CpuEulerDomainBackend::CreateState(
    const EulerDomainProblem & problem,
    const EulerDomainStateKey & key ) const
{
    ValidateEulerDomainProblem( problem );
    ValidateKey( key );

    auto state = std::make_unique< CpuEulerDomainState >();
    state->problem = problem;
    state->key = key;
    state->values.resize( ValueCount( problem ) );
    return state;
}

void CpuEulerDomainBackend::Upload(
    EulerDomainState & state,
    const EulerDomainConstFieldView & field ) const
{
    CpuEulerDomainState & cpuState = AsCpuState( state );
    ValidateEulerDomainField( cpuState.problem, field );
    std::copy(
        field.values,
        field.values + cpuState.values.size(),
        cpuState.values.begin() );
}

void CpuEulerDomainBackend::Advance(
    EulerDomainState & state,
    int steps,
    const EulerDomainRunOptions & ) const
{
    AsCpuState( state );
    if ( steps < 0 )
    {
        throw std::invalid_argument(
            "CPU Euler domain advance steps cannot be negative" );
    }
    if ( steps > 0 )
    {
        throw std::logic_error(
            "CPU Euler domain state advancement is reserved for E5" );
    }
}

void CpuEulerDomainBackend::Download(
    const EulerDomainState & state,
    EulerDomainFieldView & field ) const
{
    const CpuEulerDomainState & cpuState = AsCpuState( state );
    ValidateEulerDomainField( cpuState.problem, field );
    std::copy(
        cpuState.values.begin(),
        cpuState.values.end(),
        field.values );
}

EndNameSpace
