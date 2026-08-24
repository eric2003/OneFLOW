/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\*---------------------------------------------------------------------------*/

#include "AccelRuntime.h"
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#ifndef ONEFLOW_DEFAULT_ACCEL_BACKEND
#define ONEFLOW_DEFAULT_ACCEL_BACKEND "CPU"
#endif

BeginNameSpace( ONEFLOW )

namespace
{

bool ParseEnvironmentInteger( const char * name, int & value )
{
    const char * text = std::getenv( name );
    if ( text == nullptr || *text == '\0' ) return false;

    char * end = nullptr;
    errno = 0;
    long parsed = std::strtol( text, & end, 10 );
    if ( errno != 0 || end == text || *end != '\0'
         || parsed < 0 || parsed > INT_MAX )
    {
        return false;
    }

    value = static_cast< int >( parsed );
    return true;
}

int DetectLocalRank()
{
    const char * names[] = {
        "OMPI_COMM_WORLD_LOCAL_RANK",
        "MV2_COMM_WORLD_LOCAL_RANK",
        "SLURM_LOCALID"
    };

    int value = 0;
    for ( const char * name : names )
    {
        if ( ParseEnvironmentInteger( name, value ) ) return value;
    }
    return 0;
}

int DetectLocalSize()
{
    const char * names[] = {
        "OMPI_COMM_WORLD_LOCAL_SIZE",
        "MV2_COMM_WORLD_LOCAL_SIZE",
        "SLURM_NTASKS_PER_NODE"
    };

    int value = 1;
    for ( const char * name : names )
    {
        if ( ParseEnvironmentInteger( name, value ) ) return value;
    }
    return 1;
}

int DetectRequestedDevice()
{
    int value = -1;
    ParseEnvironmentInteger( "ONEFLOW_DEVICE_ID", value );
    return value;
}

std::string RequestedBackendName()
{
    const char * environment = std::getenv( "ONEFLOW_ACCEL_BACKEND" );
    if ( environment != nullptr && *environment != '\0' ) return environment;
    return ONEFLOW_DEFAULT_ACCEL_BACKEND;
}

}

AccelRuntime & AccelRuntime::Instance()
{
    static AccelRuntime runtime;
    return runtime;
}

void AccelRuntime::Initialize( int worldRank, int worldSize )
{
    if ( initialized ) return;

    context.worldRank = worldRank;
    context.worldSize = worldSize;
    context.localRank = DetectLocalRank();
    context.localSize = DetectLocalSize();
    context.requestedDevice = DetectRequestedDevice();
#ifdef ONEFLOW_ENABLE_MULTI_DEVICE
    context.multiDeviceEnabled = true;
#else
    context.multiDeviceEnabled = false;
#endif

    const AccelBackendKind requestedKind =
        ParseAccelBackendKind( RequestedBackendName() );

    AccelBackendRegistry & registry = AccelBackendRegistry::Instance();
    if ( ! registry.Contains( requestedKind ) )
    {
        throw std::runtime_error(
            std::string( "OneFLOW backend '" )
            + AccelBackendKindName( requestedKind )
            + "' is reserved but is not built in this configuration." );
    }

    backend = registry.Create( requestedKind );
    backend->Initialize( context );
    initialized = true;

    if ( context.worldRank == 0 )
    {
        std::cout << "OneFLOW compute backend: " << backend->Name() << "\n";
        std::cout << "OneFLOW multi-device context: "
                  << ( context.multiDeviceEnabled ? "enabled" : "disabled" )
                  << "\n";
    }
}

void AccelRuntime::Finalize()
{
    if ( ! initialized ) return;
    backend->Finalize();
    backend.reset();
    initialized = false;
}

void AccelRuntime::Synchronize()
{
    if ( initialized ) backend->Synchronize();
}

bool AccelRuntime::IsInitialized() const
{
    return initialized;
}

bool AccelRuntime::IsAccelerator() const
{
    return initialized && backend->IsAccelerator();
}

AccelBackend & AccelRuntime::Backend()
{
    if ( ! initialized )
    {
        throw std::runtime_error( "OneFLOW accelerator runtime is not initialized." );
    }
    return *backend;
}

const AccelContext & AccelRuntime::Context() const
{
    return context;
}

void InitializeAccelRuntime( int worldRank, int worldSize )
{
    AccelRuntime::Instance().Initialize( worldRank, worldSize );
}

void FinalizeAccelRuntime()
{
    AccelRuntime::Instance().Finalize();
}

EndNameSpace
