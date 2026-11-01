#include <gtest/gtest.h>

#include "AccelBackend.h"
#include "AccelRuntime.h"

#include <cstdlib>
#include <memory>
#include <stdexcept>
#include <string>

using namespace ONEFLOW;

namespace
{

int failedBackendDestructions = 0;

class FailingHipBackend final : public AccelBackend
{
public:
    ~FailingHipBackend() override
    {
        ++failedBackendDestructions;
    }

    AccelBackendKind Kind() const override
    {
        return AccelBackendKind::HIP;
    }

    const char * Name() const override
    {
        return "FAILING_HIP";
    }

    bool IsAccelerator() const override
    {
        return true;
    }

    bool SupportsMultiDevice() const override
    {
        return true;
    }

    void Initialize( AccelContext & context ) override
    {
        context.deviceCount = 8;
        context.selectedDevice = 6;
        throw std::runtime_error( "simulated backend initialization failure" );
    }

    void Finalize() override
    {
        ADD_FAILURE() << "Finalize must not be called after failed initialization";
    }

    void Synchronize() override {}
    void * Allocate( std::size_t ) override { return nullptr; }
    void Deallocate( void * ) override {}
    void Copy( void *, const void *, std::size_t, AccelCopyKind ) override {}
};

class TestCpuBackend final : public AccelBackend
{
public:
    AccelBackendKind Kind() const override
    {
        return AccelBackendKind::CPU;
    }

    const char * Name() const override
    {
        return "TEST_CPU";
    }

    bool IsAccelerator() const override
    {
        return false;
    }

    bool SupportsMultiDevice() const override
    {
        return false;
    }

    void Initialize( AccelContext & context ) override
    {
        context.deviceCount = 0;
        context.selectedDevice = -1;
    }

    void Finalize() override {}
    void Synchronize() override {}
    void * Allocate( std::size_t ) override { return nullptr; }
    void Deallocate( void * ) override {}
    void Copy( void *, const void *, std::size_t, AccelCopyKind ) override {}
};

void SetEnvironmentVariable( const char * name, const std::string & value )
{
#ifdef _WIN32
    _putenv_s( name, value.c_str() );
#else
    setenv( name, value.c_str(), 1 );
#endif
}

void ClearEnvironmentVariable( const char * name )
{
#ifdef _WIN32
    _putenv_s( name, "" );
#else
    unsetenv( name );
#endif
}

class EnvironmentVariableGuard
{
public:
    explicit EnvironmentVariableGuard( const char * name )
        : name_( name )
    {
        const char * value = std::getenv( name );
        hadValue_ = value != nullptr;
        if ( hadValue_ ) oldValue_ = value;
    }

    ~EnvironmentVariableGuard()
    {
        if ( hadValue_ )
        {
            SetEnvironmentVariable( name_.c_str(), oldValue_ );
        }
        else
        {
            ClearEnvironmentVariable( name_.c_str() );
        }
    }

private:
    std::string name_;
    std::string oldValue_;
    bool hadValue_ = false;
};

}

TEST( AccelRuntimeTest, FailedInitializationDoesNotPublishPartialState )
{
    EnvironmentVariableGuard backendSetting( "ONEFLOW_ACCEL_BACKEND" );
    SetEnvironmentVariable( "ONEFLOW_ACCEL_BACKEND", "HIP" );

    AccelBackendRegistry & registry = AccelBackendRegistry::Instance();
    registry.Register(
        AccelBackendKind::HIP,
        []() { return std::make_unique< FailingHipBackend >(); } );
    registry.Register(
        AccelBackendKind::CPU,
        []() { return std::make_unique< TestCpuBackend >(); } );

    AccelRuntime & runtime = AccelRuntime::Instance();
    ASSERT_FALSE( runtime.IsInitialized() );

    failedBackendDestructions = 0;
    EXPECT_THROW(
        runtime.Initialize( 3, 8 ),
        std::runtime_error );

    EXPECT_FALSE( runtime.IsInitialized() );
    EXPECT_EQ( runtime.Context().worldRank, 0 );
    EXPECT_EQ( runtime.Context().worldSize, 1 );
    EXPECT_EQ( runtime.Context().deviceCount, 0 );
    EXPECT_EQ( runtime.Context().selectedDevice, -1 );
    EXPECT_EQ( failedBackendDestructions, 1 );

    SetEnvironmentVariable( "ONEFLOW_ACCEL_BACKEND", "AUTO" );
    EXPECT_NO_THROW( runtime.Initialize( 3, 8 ) );

    ASSERT_TRUE( runtime.IsInitialized() );
    EXPECT_STREQ( runtime.Backend().Name(), "TEST_CPU" );
    EXPECT_EQ( runtime.Context().worldRank, 3 );
    EXPECT_EQ( runtime.Context().worldSize, 8 );
    EXPECT_EQ( runtime.Context().deviceCount, 0 );
    EXPECT_EQ( runtime.Context().selectedDevice, -1 );
    EXPECT_EQ( failedBackendDestructions, 2 );

    EXPECT_NO_THROW( runtime.Finalize() );
    EXPECT_FALSE( runtime.IsInitialized() );
}
