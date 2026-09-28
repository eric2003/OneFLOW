/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "AccelBackend.h"

#ifdef ONEFLOW_ENABLE_HIP

#include "HipKernel.h"

#include <hip/hip_runtime.h>

#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{

void CheckHip( hipError_t status, const char * operation )
{
    if ( status == hipSuccess ) return;
    throw std::runtime_error(
        std::string( operation ) + ": " + hipGetErrorString( status ) );
}

bool IsSelfTestEnabled()
{
    const char * value = std::getenv( "ONEFLOW_HIP_SELF_TEST" );
    return value != nullptr && std::strcmp( value, "1" ) == 0;
}

class HipBackend final : public AccelBackend
{
public:
    AccelBackendKind Kind() const override
    {
        return AccelBackendKind::HIP;
    }

    const char * Name() const override
    {
        return "HIP/DCU";
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
        int deviceCount = 0;
        CheckHip( hipGetDeviceCount( &deviceCount ), "hipGetDeviceCount" );
        if ( deviceCount <= 0 )
        {
            throw std::runtime_error( "HIP backend found no visible DCU device." );
        }

        context.deviceCount = deviceCount;
        const int requested = context.requestedDevice;
        const int selected =
            requested >= 0 ? requested : context.localRank % deviceCount;
        if ( selected < 0 || selected >= deviceCount )
        {
            throw std::runtime_error(
                "ONEFLOW_DEVICE_ID is outside the visible HIP device range." );
        }

        CheckHip( hipSetDevice( selected ), "hipSetDevice" );
        context.selectedDevice = selected;

        hipDeviceProp_t properties {};
        CheckHip(
            hipGetDeviceProperties( &properties, selected ),
            "hipGetDeviceProperties" );
        if ( context.worldRank == 0 )
        {
            std::printf(
                "OneFLOW HIP device: %d/%d, name=%s, gcnArch=%s\n",
                selected,
                deviceCount,
                properties.name,
                properties.gcnArchName );
        }

        if ( IsSelfTestEnabled() ) RunHipBackendSelfTest();
    }

    void Finalize() override
    {
        CheckHip( hipDeviceSynchronize(), "hipDeviceSynchronize(finalize)" );
    }

    void Synchronize() override
    {
        CheckHip( hipDeviceSynchronize(), "hipDeviceSynchronize" );
    }

    void * Allocate( std::size_t bytes ) override
    {
        if ( bytes == 0 ) return nullptr;
        void * pointer = nullptr;
        CheckHip( hipMalloc( &pointer, bytes ), "hipMalloc" );
        return pointer;
    }

    void Deallocate( void * pointer ) override
    {
        if ( pointer != nullptr ) CheckHip( hipFree( pointer ), "hipFree" );
    }

    void Copy(
        void * destination,
        const void * source,
        std::size_t bytes,
        AccelCopyKind kind ) override
    {
        if ( bytes == 0 ) return;
        hipMemcpyKind hipKind = hipMemcpyDefault;
        switch ( kind )
        {
            case AccelCopyKind::HostToHost:
                hipKind = hipMemcpyHostToHost;
                break;
            case AccelCopyKind::HostToDevice:
                hipKind = hipMemcpyHostToDevice;
                break;
            case AccelCopyKind::DeviceToHost:
                hipKind = hipMemcpyDeviceToHost;
                break;
            case AccelCopyKind::DeviceToDevice:
                hipKind = hipMemcpyDeviceToDevice;
                break;
        }
        CheckHip( hipMemcpy( destination, source, bytes, hipKind ), "hipMemcpy" );
    }
};

struct RegisterHipBackend
{
    RegisterHipBackend()
    {
        AccelBackendRegistry::Instance().Register(
            AccelBackendKind::HIP,
            []() {
                return std::unique_ptr< AccelBackend >( new HipBackend() );
            } );
    }
};

RegisterHipBackend registerHipBackend;

}

EndNameSpace

#endif
