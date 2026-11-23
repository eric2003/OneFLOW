/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "MpiTrace.h"

#include "NamespaceMacros.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <string>
#include <unordered_map>

BeginNameSpace( ONEFLOW )

namespace
{

struct Counter
{
    unsigned long long calls = 0;
    unsigned long long bytes = 0;
    double seconds = 0.0;
};

std::mutex & TraceMutex()
{
    static std::mutex mutex;
    return mutex;
}

std::unordered_map< std::string, Counter > & Counters()
{
    static std::unordered_map< std::string, Counter > counters;
    return counters;
}

double Now()
{
    using clock = std::chrono::steady_clock;
    return std::chrono::duration< double >(
        clock::now().time_since_epoch() ).count();
}

// "%r" expands to the MPI rank so ranks do not share one file.
std::string RankAwarePath( const char * configured )
{
    std::string path( configured == nullptr ? "" : configured );
    const std::string token( "%r" );
    const std::size_t position = path.find( token );
    if ( position == std::string::npos ) return path;

    const char * rankNames[] = {
        "OMPI_COMM_WORLD_RANK", "MV2_COMM_WORLD_RANK", "SLURM_PROCID" };
    const char * rank = "0";
    for ( const char * name : rankNames )
    {
        const char * value = std::getenv( name );
        if ( value != nullptr && value[ 0 ] != '\0' )
        {
            rank = value;
            break;
        }
    }
    path.replace( position, token.size(), rank );
    return path;
}

} // namespace

bool MpiTrace::Enabled()
{
    static const bool enabled = []()
    {
        const char * path = std::getenv( "ONEFLOW_MPI_TRACE_FILE" );
        return path != nullptr && path[ 0 ] != '\0';
    }();
    return enabled;
}

MpiTrace::Scope::Scope( const char * kind, std::size_t bytes )
    : kind_( kind ), bytes_( bytes ), active_( MpiTrace::Enabled() )
{
    if ( active_ ) start_ = Now();
}

MpiTrace::Scope::~Scope()
{
    if ( ! active_ ) return;
    const double elapsed = Now() - start_;
    try
    {
        std::lock_guard< std::mutex > lock( TraceMutex() );
        Counter & counter = Counters()[ std::string( kind_ ) ];
        ++ counter.calls;
        counter.bytes += bytes_;
        counter.seconds += elapsed;
    }
    catch ( ... )
    {
        // Tracing must never change solver behaviour.
    }
}

void MpiTrace::Flush()
{
    const char * configured = std::getenv( "ONEFLOW_MPI_TRACE_FILE" );
    if ( configured == nullptr || configured[ 0 ] == '\0' ) return;

    const std::string path = RankAwarePath( configured );
    if ( path.empty() ) return;

    std::lock_guard< std::mutex > lock( TraceMutex() );
    std::ofstream output( path, std::ios::out | std::ios::trunc );
    if ( ! output ) return;

    output << "kind\tcalls\tbytes\tseconds\n";
    for ( const auto & entry : Counters() )
    {
        output << entry.first << '\t' << entry.second.calls << '\t'
               << entry.second.bytes << '\t' << entry.second.seconds << '\n';
    }
}

EndNameSpace
