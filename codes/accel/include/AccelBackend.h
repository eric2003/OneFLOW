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

#pragma once

#include "NamespaceMacros.h"
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

enum class AccelBackendKind
{
    CPU,
    HIP,
    CUDA,
    KOKKOS
};

enum class AccelCopyKind
{
    HostToHost,
    HostToDevice,
    DeviceToHost,
    DeviceToDevice
};

struct AccelContext
{
    int worldRank = 0;
    int worldSize = 1;
    int localRank = 0;
    int localSize = 1;
    int requestedDevice = -1;
    int selectedDevice = -1;
    int deviceCount = 0;
    bool multiDeviceEnabled = false;
};

class AccelBackend
{
public:
    virtual ~AccelBackend() = default;

    virtual AccelBackendKind Kind() const = 0;
    virtual const char * Name() const = 0;
    virtual bool IsAccelerator() const = 0;
    virtual bool SupportsMultiDevice() const = 0;

    virtual void Initialize( AccelContext & context ) = 0;
    virtual void Finalize() = 0;
    virtual void Synchronize() = 0;

    virtual void * Allocate( std::size_t bytes ) = 0;
    virtual void Deallocate( void * pointer ) = 0;
    virtual void Copy(
        void * destination,
        const void * source,
        std::size_t bytes,
        AccelCopyKind kind ) = 0;
};

using AccelBackendFactory = std::function< std::unique_ptr< AccelBackend >() >;

class AccelBackendRegistry
{
public:
    static AccelBackendRegistry & Instance();

    void Register( AccelBackendKind kind, AccelBackendFactory factory );
    bool Contains( AccelBackendKind kind ) const;
    std::unique_ptr< AccelBackend > Create( AccelBackendKind kind ) const;
    std::vector< AccelBackendKind > RegisteredBackends() const;

private:
    AccelBackendRegistry() = default;
    AccelBackendRegistry( const AccelBackendRegistry & ) = delete;
    AccelBackendRegistry & operator=( const AccelBackendRegistry & ) = delete;

private:
    struct Entry;
    std::vector< std::unique_ptr< Entry > > entries;
};

AccelBackendKind ParseAccelBackendKind( const std::string & name );
const char * AccelBackendKindName( AccelBackendKind kind );

EndNameSpace
