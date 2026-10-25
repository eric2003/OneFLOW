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

    // Grow-only.  Several call sites request slightly different "total cell"
    // counts for the same resident buffer (the boundary model uses
    // nCells + nGhostCells, the interface upload uses nCells + nBoundaryFaces).
    // A shrink-and-reallocate made hipMalloc/hipFree run four times per step;
    // that cost 11-14 ms/step on the 8-zone DCU case, independent of the grid
    // size.  Keeping the larger allocation removes the churn and every caller
    // already guards with "Size() < needed".
    void Resize( std::size_t count )
    {
        if ( count <= count_ ) return;

        T * grown = static_cast< T * >(
            AccelRuntime::Instance().Backend().Allocate(
                count * sizeof( T ) ) );
        if ( pointer_ != nullptr )
        {
            AccelRuntime::Instance().Backend().Copy(
                grown, pointer_, count_ * sizeof( T ),
                AccelCopyKind::DeviceToDevice );
            if ( AccelRuntime::Instance().IsInitialized() )
            {
                AccelRuntime::Instance().Backend().Deallocate( pointer_ );
            }
        }
        pointer_ = grown;
        count_ = count;
    }

    void Reset()
    {
        if ( pointer_ != nullptr )
        {
            // A failure path may reach static destruction after the
            // accelerator runtime has already been finalized.
            if ( AccelRuntime::Instance().IsInitialized() )
            {
                AccelRuntime::Instance().Backend().Deallocate( pointer_ );
            }
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
