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

#include "AccelRuntime.h"
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

BeginNameSpace( ONEFLOW )

namespace DeviceBufferDetail
{

class AllocationGuard
{
public:
    AllocationGuard( AccelBackend & backend, void * pointer ) noexcept
        : backend_( backend ), pointer_( pointer )
    {
    }

    ~AllocationGuard() noexcept
    {
        if ( pointer_ == nullptr ) return;
        try
        {
            backend_.Deallocate( pointer_ );
        }
        catch ( ... )
        {
            // Destruction is best-effort. Explicit Reset() reports release errors.
        }
    }

    AllocationGuard( const AllocationGuard & ) = delete;
    AllocationGuard & operator=( const AllocationGuard & ) = delete;

    void Release() noexcept
    {
        pointer_ = nullptr;
    }

private:
    AccelBackend & backend_;
    void * pointer_ = nullptr;
};

}

template< typename T >
class DeviceBuffer
{
public:
    DeviceBuffer() = default;

    explicit DeviceBuffer( std::size_t count )
    {
        Resize( count );
    }

    ~DeviceBuffer() noexcept
    {
        ResetNoThrow();
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
            ResetNoThrow();
            Swap( other );
        }
        return *this;
    }

    // Grow-only. Different solver stages can request different total-cell
    // extents for the same resident buffer. Retaining the largest allocation
    // avoids repeated device allocation/free operations in the time loop.
    //
    // Allocation and copy failures leave the current buffer unchanged. Once
    // the replacement has been committed, failure to release the old
    // allocation is reported while the buffer continues to own the new one.
    void Resize( std::size_t count )
    {
        if ( count <= count_ ) return;
        if ( count > ( std::numeric_limits< std::size_t >::max )() / sizeof( T ) )
        {
            throw std::length_error( "DeviceBuffer size exceeds addressable memory." );
        }

        AccelRuntime & runtime = AccelRuntime::Instance();
        AccelBackend & backend = runtime.Backend();
        T * grown = static_cast< T * >(
            backend.Allocate( count * sizeof( T ) ) );
        DeviceBufferDetail::AllocationGuard allocationGuard( backend, grown );

        if ( pointer_ != nullptr )
        {
            backend.Copy(
                grown, pointer_, count_ * sizeof( T ),
                AccelCopyKind::DeviceToDevice );
        }

        T * previous = pointer_;
        pointer_ = grown;
        count_ = count;
        allocationGuard.Release();

        if ( previous != nullptr ) backend.Deallocate( previous );
    }

    void Reset()
    {
        T * previous = pointer_;
        pointer_ = nullptr;
        count_ = 0;

        if ( previous != nullptr )
        {
            AccelRuntime & runtime = AccelRuntime::Instance();
            if ( runtime.IsInitialized() )
            {
                runtime.Backend().Deallocate( previous );
            }
        }
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
    void ResetNoThrow() noexcept
    {
        try
        {
            Reset();
        }
        catch ( ... )
        {
            // A destructor or noexcept move assignment must not propagate
            // backend release failures. Reset() remains available to callers
            // that need to observe such failures.
        }
    }

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
