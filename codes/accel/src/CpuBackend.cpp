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

#include "AccelBackend.h"
#include <cstring>
#include <memory>
#include <new>

BeginNameSpace( ONEFLOW )

namespace
{

class CpuBackend final : public AccelBackend
{
public:
    AccelBackendKind Kind() const override
    {
        return AccelBackendKind::CPU;
    }

    const char * Name() const override
    {
        return "CPU";
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

    void Finalize() override
    {
    }

    void Synchronize() override
    {
    }

    void * Allocate( std::size_t bytes ) override
    {
        if ( bytes == 0 ) return nullptr;
        return ::operator new( bytes );
    }

    void Deallocate( void * pointer ) override
    {
        ::operator delete( pointer );
    }

    void Copy(
        void * destination,
        const void * source,
        std::size_t bytes,
        AccelCopyKind ) override
    {
        if ( bytes != 0 ) std::memcpy( destination, source, bytes );
    }
};

struct RegisterCpuBackend
{
    RegisterCpuBackend()
    {
        AccelBackendRegistry::Instance().Register(
            AccelBackendKind::CPU,
            []() {
                return std::unique_ptr< AccelBackend >( new CpuBackend() );
            } );
    }
};

RegisterCpuBackend registerCpuBackend;

}

EndNameSpace
