#include "OneDEulerBackend.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace oneflow_1d
{
namespace
{

void CheckProblem( const EulerProblem & problem )
{
    if ( problem.nx <= 0 || problem.gamma <= 1.0 || problem.dt <= 0.0
         || problem.dx <= 0.0 )
    {
        throw std::invalid_argument( "invalid Euler problem" );
    }
}

struct CpuEulerState final : EulerState
{
    explicit CpuEulerState( const EulerProblem & value )
        : problem( value ), values( EulerComponents * value.nx, 0.0 )
    {
    }

    EulerProblem problem;
    std::vector< double > values;
    bool uploaded = false;
};

CpuEulerState & CpuState( EulerState & state )
{
    auto * result = dynamic_cast< CpuEulerState * >( &state );
    if ( result == nullptr ) throw std::invalid_argument( "state is not a CPU Euler state" );
    return *result;
}

const CpuEulerState & CpuState( const EulerState & state )
{
    auto * result = dynamic_cast< const CpuEulerState * >( &state );
    if ( result == nullptr ) throw std::invalid_argument( "state is not a CPU Euler state" );
    return *result;
}

void CheckSteps( int steps, const EulerRunOptions & options )
{
    if ( steps <= 0 ) throw std::invalid_argument( "Euler steps must be positive" );
    if ( options.mode == EulerRunMode::FullTrace
         && ( steps != 1 || options.trace == nullptr ) )
    {
        throw std::invalid_argument(
            "FullTrace requires exactly one step and a trace" );
    }
    if ( options.mode == EulerRunMode::NoTrace && options.trace != nullptr )
    {
        throw std::invalid_argument( "NoTrace cannot provide a trace" );
    }
}

} // namespace

void EulerBackend::Step(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace ) const
{
    EulerProblem problem{ nx, gamma, dt, dx, boundary };
    auto execution = CreateState( problem );
    Upload( *execution, state );
    Advance( *execution, 1, { EulerRunMode::FullTrace, &trace, nullptr } );
}

const char * CpuEulerBackend::Name() const
{
    return "CPU";
}

bool CpuEulerBackend::IsAccelerator() const
{
    return false;
}

std::unique_ptr< EulerState > CpuEulerBackend::CreateState(
    const EulerProblem & problem ) const
{
    CheckProblem( problem );
    return std::make_unique< CpuEulerState >( problem );
}

void CpuEulerBackend::Upload(
    EulerState & state, const double * hostState ) const
{
    auto & cpu = CpuState( state );
    if ( hostState == nullptr ) throw std::invalid_argument( "null Euler host state" );
    std::copy( hostState, hostState + cpu.values.size(), cpu.values.begin() );
    cpu.uploaded = true;
}

void CpuEulerBackend::Advance(
    EulerState & state, int steps, const EulerRunOptions & options ) const
{
    CheckSteps( steps, options );
    auto & cpu = CpuState( state );
    if ( ! cpu.uploaded ) throw std::logic_error( "CPU Euler state was not uploaded" );
    EulerTrace localTrace;
    EulerTrace * trace = options.mode == EulerRunMode::FullTrace
        ? options.trace : &localTrace;
    for ( int step = 0; step < steps; ++ step )
    {
        OneDCpuEulerStep(
            cpu.values.data(), cpu.problem.nx, cpu.problem.gamma,
            cpu.problem.dt, cpu.problem.dx, cpu.problem.boundary, *trace );
        std::copy(
            trace->state.begin() + EulerRkStages * cpu.values.size(),
            trace->state.end(), cpu.values.begin() );
    }
}

void CpuEulerBackend::Download(
    const EulerState & state, double * hostState ) const
{
    const auto & cpu = CpuState( state );
    if ( ! cpu.uploaded ) throw std::logic_error( "CPU Euler state was not uploaded" );
    if ( hostState == nullptr ) throw std::invalid_argument( "null Euler host state" );
    std::copy( cpu.values.begin(), cpu.values.end(), hostState );
}

void CpuEulerBackend::Step(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace ) const
{
    EulerBackend::Step( state, nx, gamma, dt, dx, boundary, trace );
}

} // namespace oneflow_1d
