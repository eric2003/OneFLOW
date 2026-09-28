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
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>

BeginNameSpace( ONEFLOW )

struct AccelBackendRegistry::Entry
{
    AccelBackendKind kind;
    AccelBackendFactory factory;
};

AccelBackendRegistry & AccelBackendRegistry::Instance()
{
    static AccelBackendRegistry registry;
    return registry;
}

void AccelBackendRegistry::Register(
    AccelBackendKind kind,
    AccelBackendFactory factory )
{
    for ( auto & entry : entries )
    {
        if ( entry->kind == kind )
        {
            entry->factory = std::move( factory );
            return;
        }
    }

    std::unique_ptr< Entry > entry( new Entry() );
    entry->kind = kind;
    entry->factory = std::move( factory );
    entries.push_back( std::move( entry ) );
}

bool AccelBackendRegistry::Contains( AccelBackendKind kind ) const
{
    for ( const auto & entry : entries )
    {
        if ( entry->kind == kind ) return true;
    }
    return false;
}

std::unique_ptr< AccelBackend >
AccelBackendRegistry::Create( AccelBackendKind kind ) const
{
    for ( const auto & entry : entries )
    {
        if ( entry->kind == kind ) return entry->factory();
    }

    throw std::runtime_error(
        std::string( "OneFLOW backend is not registered: " )
        + AccelBackendKindName( kind ) );
}

std::vector< AccelBackendKind >
AccelBackendRegistry::RegisteredBackends() const
{
    std::vector< AccelBackendKind > result;
    result.reserve( entries.size() );
    for ( const auto & entry : entries ) result.push_back( entry->kind );
    return result;
}

AccelBackendKind ParseAccelBackendKind( const std::string & name )
{
    std::string normalized = name;
    std::transform(
        normalized.begin(), normalized.end(), normalized.begin(),
        []( unsigned char value ) {
            return static_cast< char >( std::toupper( value ) );
        } );

    if ( normalized == "CPU" ) return AccelBackendKind::CPU;
    if ( normalized == "HIP" ) return AccelBackendKind::HIP;
    if ( normalized == "CUDA" ) return AccelBackendKind::CUDA;
    if ( normalized == "KOKKOS" ) return AccelBackendKind::KOKKOS;

    throw std::invalid_argument(
        "Unknown OneFLOW backend '" + name
        + "'. Expected CPU, HIP, CUDA, or KOKKOS." );
}

const char * AccelBackendKindName( AccelBackendKind kind )
{
    switch ( kind )
    {
        case AccelBackendKind::CPU: return "CPU";
        case AccelBackendKind::HIP: return "HIP";
        case AccelBackendKind::CUDA: return "CUDA";
        case AccelBackendKind::KOKKOS: return "KOKKOS";
    }
    return "UNKNOWN";
}

EndNameSpace
