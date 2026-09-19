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

#pragma once

#include "EulerDomain.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <unordered_map>

BeginNameSpace( ONEFLOW )

struct EulerDomainStateKeyHash
{
    std::size_t operator()( const EulerDomainStateKey & key ) const noexcept;
};

// Owns opaque backend state without exposing MRField, Zone, or backend
// implementation details to the solver. The owner of the registry decides
// when a state is created (after initialization) and when it is released
// (before finalization or solver teardown).
class EulerDomainStateRegistry
{
public:
    EulerDomainStateRegistry() = default;
    ~EulerDomainStateRegistry() = default;

    EulerDomainStateRegistry( const EulerDomainStateRegistry & ) = delete;
    EulerDomainStateRegistry( EulerDomainStateRegistry && ) = default;
    EulerDomainStateRegistry & operator=( EulerDomainStateRegistry && ) = default;
    EulerDomainStateRegistry & operator=( const EulerDomainStateRegistry & ) = delete;

    using StateFactory = std::function< std::unique_ptr< EulerDomainState >() >;

    bool Contains( const EulerDomainStateKey & key ) const;
    EulerDomainState & GetOrCreate(
        const EulerDomainStateKey & key, const StateFactory & factory );
    bool Invalidate( const EulerDomainStateKey & key );
    void Insert(
        const EulerDomainStateKey & key,
        std::unique_ptr< EulerDomainState > state );
    EulerDomainState & Get( const EulerDomainStateKey & key );
    const EulerDomainState & Get( const EulerDomainStateKey & key ) const;
    void Erase( const EulerDomainStateKey & key );
    void Clear();
    std::size_t Size() const;

private:
    using StateMap = std::unordered_map<
        EulerDomainStateKey, std::unique_ptr< EulerDomainState >,
        EulerDomainStateKeyHash >;
    StateMap states;
};

EndNameSpace
