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

#include "AccelBackend.h"
#include <memory>

BeginNameSpace( ONEFLOW )

class AccelRuntime
{
public:
    static AccelRuntime & Instance();

    void Initialize( int worldRank, int worldSize );
    void Finalize();
    void Synchronize();

    bool IsInitialized() const;
    bool IsAccelerator() const;
    AccelBackend & Backend();
    const AccelContext & Context() const;

private:
    AccelRuntime() = default;
    AccelRuntime( const AccelRuntime & ) = delete;
    AccelRuntime & operator=( const AccelRuntime & ) = delete;

private:
    AccelContext context;
    std::unique_ptr< AccelBackend > backend;
    bool initialized = false;
};

void InitializeAccelRuntime( int worldRank, int worldSize );
void FinalizeAccelRuntime();

EndNameSpace
