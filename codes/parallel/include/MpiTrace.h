/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#pragma once

#include <cstddef>

#include "NamespaceMacros.h"

BeginNameSpace( ONEFLOW )

// Opt-in accounting of every MPI call the solver makes: count, bytes and wall
// time per operation kind, flushed to "<path>.<rank>.tsv" at HXFinalize.
//
// Scaling questions ("is the halo latency- or bandwidth-bound?", "how many
// collectives run per step?") cannot be answered from stage timers alone: one
// category merges every message of every neighbour, and the deferred drain
// hides the wait inside whatever blocks next.  Counting the calls separates
// message count (latency) from bytes (bandwidth) and names the collectives.
//
// Disabled unless ONEFLOW_MPI_TRACE_FILE is set, in which case the cost is one
// cached environment lookup per call.
class MpiTrace
{
public:
    static bool Enabled();

    // RAII timer: records one call of `kind` transferring `bytes`.
    class Scope
    {
    public:
        Scope( const char * kind, std::size_t bytes );
        ~Scope();

        Scope( const Scope & ) = delete;
        Scope & operator=( const Scope & ) = delete;

    private:
        const char * kind_ = nullptr;
        std::size_t bytes_ = 0;
        bool active_ = false;
        double start_ = 0.0;
    };

    static void Flush();
};

EndNameSpace
