#include "EulerDomainStateLifecycle.h"

#include <gtest/gtest.h>

#include <stdexcept>
#include <vector>

namespace
{

using namespace ONEFLOW;

struct LifecycleState final : EulerDomainState
{
    int generation = 0;
    std::vector< Real > uploaded;
};

class MockBackend final : public EulerDomainBackend
{
public:
    explicit MockBackend( bool failUpload = false )
        : failUpload_( failUpload )
    {
    }

    const char* Name() const override { return "mock"; }
    bool IsAccelerator() const override { return false; }

    std::unique_ptr< EulerDomainState > CreateState(
        const EulerDomainProblem&, const EulerDomainStateKey& ) const override
    {
        auto state = std::make_unique< LifecycleState >();
        state->generation = ++creations;
        return state;
    }

    void Upload(
        EulerDomainState& state,
        const EulerDomainConstFieldView& field ) const override
    {
        ++uploads;
        if ( failUpload_ )
        {
            throw std::runtime_error( "mock upload failed" );
        }

        auto& lifecycleState = dynamic_cast< LifecycleState& >( state );
        lifecycleState.uploaded.assign(
            field.values,
            field.values + field.nCells * field.nEquations );
    }

    void Advance(
        EulerDomainState&, int, const EulerDomainRunOptions& ) const override
    {
    }

    void Download(
        const EulerDomainState&, EulerDomainFieldView& ) const override
    {
    }

    mutable int creations = 0;
    mutable int uploads = 0;

private:
    bool failUpload_ = false;
};

EulerDomainProblem Problem()
{
    EulerDomainProblem result;
    result.nCells = 4;
    result.nEquations = 5;
    result.dt = 0.001;
    result.dx = 0.25;
    return result;
}

EulerDomainStateKey Key()
{
    return { 0, 0, 0, AccelBackendKind::CPU };
}

TEST( EulerDomainStateLifecycle, InitializeCreatesUploadsAndReplaces )
{
    EulerDomainStateRegistry registry;
    MockBackend backend;
    const EulerDomainProblem problem = Problem();
    Real values[ 20 ] = {};
    values[ 0 ] = 1.0;
    EulerDomainConstFieldView field{ 4, 5, values };

    EulerDomainState& first = EulerDomainStateLifecycle::Initialize(
        registry, backend, problem, Key(), field );
    EXPECT_EQ( static_cast< LifecycleState& >( first ).generation, 1 );
    EXPECT_EQ( backend.creations, 1 );
    EXPECT_EQ( backend.uploads, 1 );

    EulerDomainState& second = EulerDomainStateLifecycle::Initialize(
        registry, backend, problem, Key(), field );
    EXPECT_EQ( static_cast< LifecycleState& >( second ).generation, 2 );
    EXPECT_EQ( registry.Size(), 1u );
    EXPECT_EQ( backend.creations, 2 );
    EXPECT_EQ( backend.uploads, 2 );
}

TEST( EulerDomainStateLifecycle, RestartReplacesStateAfterInvalidation )
{
    EulerDomainStateRegistry registry;
    MockBackend backend;
    const EulerDomainProblem problem = Problem();
    Real values[ 20 ] = {};
    EulerDomainConstFieldView field{ 4, 5, values };

    EulerDomainStateLifecycle::Initialize(
        registry, backend, problem, Key(), field );
    EulerDomainState& restarted = EulerDomainStateLifecycle::Restart(
        registry, backend, problem, Key(), field );

    EXPECT_EQ( static_cast< LifecycleState& >( restarted ).generation, 2 );
    EXPECT_EQ( registry.Size(), 1u );
    EXPECT_EQ( backend.creations, 2 );
    EXPECT_EQ( backend.uploads, 2 );
}

TEST( EulerDomainStateLifecycle, FailedUploadDoesNotRegisterState )
{
    EulerDomainStateRegistry registry;
    MockBackend backend( true );
    Real values[ 20 ] = {};
    EulerDomainConstFieldView field{ 4, 5, values };

    EXPECT_THROW(
        EulerDomainStateLifecycle::Initialize(
            registry, backend, Problem(), Key(), field ),
        std::runtime_error );
    EXPECT_EQ( registry.Size(), 0u );
    EXPECT_EQ( backend.creations, 1 );
    EXPECT_EQ( backend.uploads, 1 );
}

} // namespace
