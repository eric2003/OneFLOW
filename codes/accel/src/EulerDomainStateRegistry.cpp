/*---------------------------------------------------------------------------*\\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
\\*---------------------------------------------------------------------------*/

#include "EulerDomainStateRegistry.h"

#include <stdexcept>
#include <type_traits>
#include <utility>

BeginNameSpace( ONEFLOW )

namespace
{

std::size_t Mix( std::size_t seed, std::size_t value )
{
    return seed ^ ( value + static_cast< std::size_t >( 0x9e3779b9 )
        + ( seed << 6 ) + ( seed >> 2 ) );
}

}

std::size_t EulerDomainStateKeyHash::operator()(
    const EulerDomainStateKey & key ) const noexcept
{
    std::size_t result = std::hash< int >{}( key.solverIndex );
    result = Mix( result, std::hash< int >{}( key.localZoneId ) );
    result = Mix( result, std::hash< int >{}( key.gridLevel ) );
    result = Mix( result, std::hash< int >{}(
        static_cast< int >( key.backend ) ) );
    result = Mix( result, std::hash< int >{}( key.deviceId ) );
    return result;
}

bool EulerDomainStateRegistry::Contains( const EulerDomainStateKey & key ) const
{
    return states.find( key ) != states.end();
}

EulerDomainState & EulerDomainStateRegistry::GetOrCreate(
    const EulerDomainStateKey & key, const StateFactory & factory )
{
    const auto iter = states.find( key );
    if ( iter != states.end() )
    {
        return *iter->second;
    }
    if ( ! factory )
    {
        throw std::invalid_argument(
            "cannot create Euler domain state without a factory" );
    }
    Insert( key, factory() );
    return Get( key );
}

bool EulerDomainStateRegistry::Invalidate( const EulerDomainStateKey & key )
{
    return states.erase( key ) != 0;
}

void EulerDomainStateRegistry::Insert(
    const EulerDomainStateKey & key,
    std::unique_ptr< EulerDomainState > state )
{
    if ( state == nullptr )
    {
        throw std::invalid_argument( "cannot register a null Euler domain state" );
    }
    const auto result = states.emplace( key, std::move( state ) );
    if ( ! result.second )
    {
        throw std::logic_error( "Euler domain state key is already registered" );
    }
}

EulerDomainState & EulerDomainStateRegistry::Get(
    const EulerDomainStateKey & key )
{
    const auto iter = states.find( key );
    if ( iter == states.end() )
    {
        throw std::out_of_range( "Euler domain state key is not registered" );
    }
    return *iter->second;
}

const EulerDomainState & EulerDomainStateRegistry::Get(
    const EulerDomainStateKey & key ) const
{
    const auto iter = states.find( key );
    if ( iter == states.end() )
    {
        throw std::out_of_range( "Euler domain state key is not registered" );
    }
    return *iter->second;
}

void EulerDomainStateRegistry::Erase( const EulerDomainStateKey & key )
{
    states.erase( key );
}

void EulerDomainStateRegistry::Clear()
{
    states.clear();
}

std::size_t EulerDomainStateRegistry::Size() const
{
    return states.size();
}

EulerDomainStateRegistry::EulerDomainStateRegistry(
    EulerDomainStateRegistry && other ) noexcept
{
    static_assert(
        std::is_nothrow_default_constructible<StateMap>::value,
        "Accelerator state map must be nothrow default constructible");
    static_assert(
        noexcept( std::declval<StateMap&>().swap(
            std::declval<StateMap&>() ) ),
        "Accelerator state map swap must be nothrow");
    states.swap( other.states );
}


EndNameSpace
