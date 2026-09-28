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

#include "AccelRuntime.h"
#include <cstddef>
#include <utility>

BeginNameSpace( ONEFLOW )

template< typename T >
class DeviceBuffer
{
public:
    DeviceBuffer() = default;

    explicit DeviceBuffer( std::size_t count )
    {
        Resize( count );
    }

    ~DeviceBuffer()
    {
        Reset();
    }

    DeviceBuffer( const DeviceBuffer & ) = delete;
    DeviceBuffer & operator=( const DeviceBuffer & ) = delete;

    DeviceBuffer( DeviceBuffer && other ) noexcept
    {
        Swap( other );
    }

    DeviceBuffer & operator=( DeviceBuffer && other ) noexcept
    {
        if ( this != & other )
        {
            Reset();
            Swap( other );
        }
        return *this;
    }

    void Resize( std::size_t count )
    {
        if ( count == count_ ) return;
        Reset();
        if ( count == 0 ) return;

        pointer_ = static_cast< T * >(
            AccelRuntime::Instance().Backend().Allocate(
                count * sizeof( T ) ) );
        count_ = count;
    }

    void Reset()
    {
        if ( pointer_ != nullptr )
        {
            AccelRuntime::Instance().Backend().Deallocate( pointer_ );
        }
        pointer_ = nullptr;
        count_ = 0;
    }

    T * Data()
    {
        return pointer_;
    }

    const T * Data() const
    {
        return pointer_;
    }

    std::size_t Size() const
    {
        return count_;
    }

    void CopyFromHost( const T * source, std::size_t count )
    {
        Resize( count );
        AccelRuntime::Instance().Backend().Copy(
            pointer_, source, count * sizeof( T ), AccelCopyKind::HostToDevice );
    }

    void CopyToHost( T * destination, std::size_t count ) const
    {
        AccelRuntime::Instance().Backend().Copy(
            destination, pointer_, count * sizeof( T ), AccelCopyKind::DeviceToHost );
    }

private:
    void Swap( DeviceBuffer & other ) noexcept
    {
        std::swap( pointer_, other.pointer_ );
        std::swap( count_, other.count_ );
    }

private:
    T * pointer_ = nullptr;
    std::size_t count_ = 0;
};

EndNameSpace
