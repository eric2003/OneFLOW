#include "CpuEulerDomainBackend.h"
#include "EulerDomainStateLifecycle.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace
{

using namespace ONEFLOW;

void RecordStage( int step, int stage, void * userData )
{
    auto * events =
        static_cast< std::vector< std::pair< int, int > > * >( userData );
    events->emplace_back( step, stage );
}

EulerDomainProblem Problem()
{
    EulerDomainProblem result;
    result.nCells = 3;
    result.nEquations = 5;
    result.dt = 0.001;
    result.dx = 0.25;
    return result;
}

EulerDomainStateKey Key()
{
    return { 2, 7, 1, AccelBackendKind::CPU };
}

TEST( CpuEulerDomainBackend, UploadsAndDownloadsEquationMajorState )
{
    CpuEulerDomainBackend backend;
    EulerDomainStateRegistry registry;
    const EulerDomainProblem problem = Problem();
    std::vector< Real > input( 15 );
    for ( int index = 0; index < 15; ++ index )
    {
        input[ index ] = static_cast< Real >( index + 1 );
    }
    EulerDomainConstFieldView source{ 3, 5, input.data() };

    EulerDomainState & state = EulerDomainStateLifecycle::Initialize(
        registry, backend, problem, Key(), source );

    std::vector< Real > output( 15, 0.0 );
    EulerDomainFieldView target{ 3, 5, output.data() };
    backend.Download( state, target );

    EXPECT_EQ( output, input );
    EXPECT_EQ( backend.Name(), std::string( "CPU-EulerDomain" ) );
    EXPECT_FALSE( backend.IsAccelerator() );
}

TEST( CpuEulerDomainBackend, RejectsNonCpuKey )
{
    CpuEulerDomainBackend backend;
    EulerDomainStateKey hipKey{ 0, 0, 0, AccelBackendKind::HIP, 0 };

    EXPECT_THROW(
        backend.CreateState( Problem(), hipKey ), std::invalid_argument );
}

TEST( CpuEulerDomainBackend, SchedulesUploadedStateStages )
{
    CpuEulerDomainBackend backend;
    EulerDomainStateRegistry registry;
    const EulerDomainProblem problem = Problem();
    Real values[ 15 ] = {};
    EulerDomainConstFieldView source{ 3, 5, values };
    EulerDomainState & state = EulerDomainStateLifecycle::Initialize(
        registry, backend, problem, Key(), source );

    std::vector< std::pair< int, int > > events;
    EulerDomainRunOptions options;
    options.stageCount = 3;
    options.stageCallback = & RecordStage;
    options.stageContext = & events;

    EXPECT_NO_THROW( backend.Advance( state, 0, options ) );
    EXPECT_NO_THROW( backend.Advance( state, 2, options ) );
    EXPECT_EQ(
        events,
        ( std::vector< std::pair< int, int > >{
            { 0, 0 }, { 0, 1 }, { 0, 2 },
            { 1, 0 }, { 1, 1 }, { 1, 2 } } ) );
}

TEST( CpuEulerDomainBackend, RequiresCallbackForPositiveAdvance )
{
    CpuEulerDomainBackend backend;
    EulerDomainStateRegistry registry;
    const EulerDomainProblem problem = Problem();
    Real values[ 15 ] = {};
    EulerDomainConstFieldView source{ 3, 5, values };
    EulerDomainState & state = EulerDomainStateLifecycle::Initialize(
        registry, backend, problem, Key(), source );

    EulerDomainRunOptions options;
    EXPECT_THROW(
        backend.Advance( state, 1, options ), std::logic_error );
}

} // namespace
