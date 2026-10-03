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

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#pragma once
#include "NamespaceMacros.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class PointerWrap
{
public:
    PointerWrap() {};
    virtual ~PointerWrap() {};
public:
    virtual void * GetPointer() { return 0; };
};

template < typename T >
class DataPointer : public PointerWrap
{
public:
    DataPointer()
        : data( std::make_unique<T>() )
    {
    }

    // Adopt ownership of an existing T*.
    explicit DataPointer( T * ptr )
        : data( ptr )
    {
    }

    // Adopt ownership from unique_ptr.
    explicit DataPointer( std::unique_ptr<T> ptr )
        : data( std::move( ptr ) )
    {
    }

    ~DataPointer() override = default;
protected:
    std::unique_ptr<T> data;
public:
    void * GetPointer() override { return data.get(); };
};

EndNameSpace
