#include <gtest/gtest.h>

#include "AccelBackend.h"
#include "AccelRuntime.h"
#include "DeviceBuffer.h"

#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <new>
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

int trackedAllocations = 0;
int trackedDeallocations = 0;
bool failDeviceCopy = false;
bool throwAfterDeviceDeallocation = false;

class TrackingCpuBackend final : public AccelBackend
{
public:
    AccelBackendKind Kind() const override { return AccelBackendKind::CPU; }
    const char * Name() const override { return "TRACKING_CPU"; }
    bool IsAccelerator() const override { return false; }
    bool SupportsMultiDevice() const override { return false; }
    void Initialize( AccelContext & context ) override
    {
        context.deviceCount = 0;
        context.selectedDevice = -1;
    }
    void Finalize() override {}
    void Synchronize() override {}
    void * Allocate( std::size_t bytes ) override
    {
        ++trackedAllocations;
        return bytes == 0 ? nullptr : ::operator new( bytes );
    }
    void Deallocate( void * pointer ) override
    {
        ++trackedDeallocations;
        ::operator delete( pointer );
        if ( throwAfterDeviceDeallocation )
        {
            throw std::runtime_error( "simulated deallocation diagnostic failure" );
        }
    }
    void Copy(
        void * destination,
        const void * source,
        std::size_t bytes,
        AccelCopyKind ) override
    {
        if ( failDeviceCopy )
        {
            throw std::runtime_error( "simulated device copy failure" );
        }
        if ( bytes != 0 ) std::memcpy( destination, source, bytes );
    }
};

class AccelRuntimeGuard
{
public:
    ~AccelRuntimeGuard() noexcept
    {
        AccelRuntime & runtime = AccelRuntime::Instance();
        if ( runtime.IsInitialized() )
        {
            try
            {
                runtime.Finalize();
            }
            catch ( ... )
            {
            }
        }
    }
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

TEST( DeviceBufferTest, FailedGrowthReleasesCandidateAndPreservesCurrentBuffer )
{
    EnvironmentVariableGuard backendSetting( "ONEFLOW_ACCEL_BACKEND" );
    SetEnvironmentVariable( "ONEFLOW_ACCEL_BACKEND", "CPU" );

    trackedAllocations = 0;
    trackedDeallocations = 0;
    failDeviceCopy = false;
    throwAfterDeviceDeallocation = false;

    AccelBackendRegistry::Instance().Register(
        AccelBackendKind::CPU,
        []() { return std::make_unique< TrackingCpuBackend >(); } );

    AccelRuntimeGuard runtimeGuard;
    AccelRuntime & runtime = AccelRuntime::Instance();
    ASSERT_FALSE( runtime.IsInitialized() );
    ASSERT_NO_THROW( runtime.Initialize( 0, 1 ) );

    DeviceBuffer< int > buffer( 2 );
    int * originalData = buffer.Data();
    ASSERT_NE( originalData, nullptr );
    ASSERT_EQ( buffer.Size(), 2u );

    failDeviceCopy = true;
    EXPECT_THROW( buffer.Resize( 4 ), std::runtime_error );
    failDeviceCopy = false;

    EXPECT_EQ( buffer.Data(), originalData );
    EXPECT_EQ( buffer.Size(), 2u );
    EXPECT_EQ( trackedAllocations, 2 );
    EXPECT_EQ( trackedDeallocations, 1 );
    EXPECT_THROW(
        buffer.Resize( ( std::numeric_limits< std::size_t >::max )() ),
        std::length_error );
    int hostValues[ 3 ] = {};
    EXPECT_THROW( buffer.CopyToHost( hostValues, 3 ), std::out_of_range );
    EXPECT_EQ( buffer.Data(), originalData );
    EXPECT_EQ( buffer.Size(), 2u );
    EXPECT_EQ( trackedAllocations, 2 );

    EXPECT_NO_THROW( buffer.Reset() );
    EXPECT_EQ( trackedDeallocations, 2 );
}

TEST( DeviceBufferTest, DestructorDoesNotPropagateBackendReleaseFailure )
{
    EnvironmentVariableGuard backendSetting( "ONEFLOW_ACCEL_BACKEND" );
    SetEnvironmentVariable( "ONEFLOW_ACCEL_BACKEND", "CPU" );

    trackedAllocations = 0;
    trackedDeallocations = 0;
    failDeviceCopy = false;
    throwAfterDeviceDeallocation = false;

    AccelBackendRegistry::Instance().Register(
        AccelBackendKind::CPU,
        []() { return std::make_unique< TrackingCpuBackend >(); } );

    AccelRuntimeGuard runtimeGuard;
    AccelRuntime & runtime = AccelRuntime::Instance();
    ASSERT_FALSE( runtime.IsInitialized() );
    ASSERT_NO_THROW( runtime.Initialize( 0, 1 ) );

    throwAfterDeviceDeallocation = true;
    EXPECT_NO_THROW(
        {
            DeviceBuffer< int > buffer( 1 );
        } );
    throwAfterDeviceDeallocation = false;

    EXPECT_EQ( trackedAllocations, 1 );
    EXPECT_EQ( trackedDeallocations, 1 );
}

TEST( DeviceBufferTest, ReleaseFailureKeepsCommittedBufferStateConsistent )
{
    EnvironmentVariableGuard backendSetting( "ONEFLOW_ACCEL_BACKEND" );
    SetEnvironmentVariable( "ONEFLOW_ACCEL_BACKEND", "CPU" );

    trackedAllocations = 0;
    trackedDeallocations = 0;
    failDeviceCopy = false;
    throwAfterDeviceDeallocation = false;

    AccelBackendRegistry::Instance().Register(
        AccelBackendKind::CPU,
        []() { return std::make_unique< TrackingCpuBackend >(); } );

    AccelRuntimeGuard runtimeGuard;
    AccelRuntime & runtime = AccelRuntime::Instance();
    ASSERT_FALSE( runtime.IsInitialized() );
    ASSERT_NO_THROW( runtime.Initialize( 0, 1 ) );

    DeviceBuffer< int > buffer( 2 );
    ASSERT_NE( buffer.Data(), nullptr );

    throwAfterDeviceDeallocation = true;
    EXPECT_THROW( buffer.Resize( 4 ), std::runtime_error );
    throwAfterDeviceDeallocation = false;

    EXPECT_NE( buffer.Data(), nullptr );
    EXPECT_EQ( buffer.Size(), 4u );
    EXPECT_EQ( trackedAllocations, 2 );
    EXPECT_EQ( trackedDeallocations, 1 );

    throwAfterDeviceDeallocation = true;
    EXPECT_THROW( buffer.Reset(), std::runtime_error );
    throwAfterDeviceDeallocation = false;

    EXPECT_EQ( buffer.Data(), nullptr );
    EXPECT_EQ( buffer.Size(), 0u );
    EXPECT_EQ( trackedDeallocations, 2 );
}

TEST( DeviceBufferTest, MoveAssignmentDoesNotPropagateBackendReleaseFailure )
{
    EnvironmentVariableGuard backendSetting( "ONEFLOW_ACCEL_BACKEND" );
    SetEnvironmentVariable( "ONEFLOW_ACCEL_BACKEND", "CPU" );

    trackedAllocations = 0;
    trackedDeallocations = 0;
    failDeviceCopy = false;
    throwAfterDeviceDeallocation = false;

    AccelBackendRegistry::Instance().Register(
        AccelBackendKind::CPU,
        []() { return std::make_unique< TrackingCpuBackend >(); } );

    AccelRuntimeGuard runtimeGuard;
    AccelRuntime & runtime = AccelRuntime::Instance();
    ASSERT_FALSE( runtime.IsInitialized() );
    ASSERT_NO_THROW( runtime.Initialize( 0, 1 ) );

    DeviceBuffer< int > destination( 1 );
    DeviceBuffer< int > source( 2 );

    throwAfterDeviceDeallocation = true;
    EXPECT_NO_THROW( destination = std::move( source ) );
    throwAfterDeviceDeallocation = false;

    EXPECT_EQ( destination.Size(), 2u );
    EXPECT_NE( destination.Data(), nullptr );
    EXPECT_EQ( source.Size(), 0u );
    EXPECT_EQ( source.Data(), nullptr );
    EXPECT_EQ( trackedDeallocations, 1 );

    EXPECT_NO_THROW( destination.Reset() );
    EXPECT_EQ( trackedDeallocations, 2 );
}
