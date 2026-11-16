/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
-------------------------------------------------------------------------------
    Opt-in stage timing for accelerator/dataflow investigations.
\*---------------------------------------------------------------------------*/

#include "StageProfiler.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

BeginNameSpace( ONEFLOW )

namespace
{

std::mutex & OutputMutex()
{
    static std::mutex mutex;
    return mutex;
}

struct CategoryTotal
{
    const char * name;
    double elapsedMs;
};

struct TaskTotal
{
    std::string name;
    unsigned long long count = 0;
    double elapsedMs = 0.0;
};

std::unordered_map< int, TaskTotal > & TaskTotals()
{
    static std::unordered_map< int, TaskTotal > totals;
    return totals;
}

struct HostTotal
{
    unsigned long long count = 0;
    double inclusiveMs = 0.0;
    double exclusiveMs = 0.0;
};

std::map< std::string, HostTotal > & HostTotals()
{
    static std::map< std::string, HostTotal > totals;
    return totals;
}

unsigned long long & DroppedHostRecords()
{
    static unsigned long long count = 0;
    return count;
}

bool HostAccountingEnabled()
{
    static const bool requested = []()
    {
        const char * value = std::getenv( "ONEFLOW_HOST_ACCOUNTING" );
        return value != nullptr && std::strcmp( value, "1" ) == 0;
    }();
    return requested && StageProfiler::Enabled();
}

auto & Totals()
{
    // Infer the size so adding a category cannot leave null trailing slots.
    static CategoryTotal totals[] = {
        { "initialization", 0.0 },
        { "gradient", 0.0 },
        { "limiter", 0.0 },
        { "reconstruction", 0.0 },
        { "boundary_reconstruction", 0.0 },
        { "host_pack", 0.0 },
        { "flux_backend", 0.0 },
        { "geometry_connectivity_H2D", 0.0 },
        { "H2D", 0.0 },
        { "HIP kernel", 0.0 },
        { "D2H", 0.0 },
        { "viscous_turbulence", 0.0 },
        { "rk_update", 0.0 },
        { "residual_state_update", 0.0 },
        { "mpi_interface", 0.0 },
        { "output", 0.0 },
        { "load_q", 0.0 },
        { "calc_time_step", 0.0 },
        { "load_residuals", 0.0 },
        { "update_residuals", 0.0 },
        { "calc_lhs", 0.0 },
        { "update_flowfield", 0.0 },
        { "calc_boundary", 0.0 },
        { "boundary_gamma_inner", 0.0 },
        { "boundary_viscosity_inner", 0.0 },
        { "boundary_bc", 0.0 },
        { "boundary_gamma_ghost", 0.0 },
        { "boundary_viscosity_ghost", 0.0 },
        { "boundary_set_id", 0.0 },
        { "boundary_prepare_data", 0.0 },
        { "boundary_calc_face_bc", 0.0 },
        { "boundary_update_bc", 0.0 },
        { "gradient_download", 0.0 },
        { "gradient_kernel", 0.0 },
        { "halo_ghost_sync", 0.0 },
        { "state_update", 0.0 },
        { "gradient_mirror", 0.0 },
        { "post_interface", 0.0 },
        { "post_residual_gate", 0.0 },
        { "post_force_gate", 0.0 },
        { "post_unsteady_dispatch", 0.0 },
        { "rhs_inviscid", 0.0 },
        { "rhs_viscous", 0.0 },
        { "rhs_source", 0.0 },
        { "flux_assembly", 0.0 },
        { "residual_metrics", 0.0 },
        { "residual_zero", 0.0 },
        { "host_mirror_refresh", 0.0 },
        { "boundary_device", 0.0 },
        { "swap_upload", 0.0 },
        { "swap_download", 0.0 },
        { "swap_store", 0.0 },
        { "task_plan", 0.0 },
        { "interface_host_send", 0.0 },
        { "interface_host_mpi", 0.0 },
        { "interface_host_recv", 0.0 },
        { "interface_ready", 0.0 },
        { "interface_upload", 0.0 },
        { "interface_pack", 0.0 },
        { "interface_mpi", 0.0 },
        { "interface_unpack", 0.0 },
        { "interface_ghost_sync", 0.0 },
        { "post_upload", 0.0 },
        { "post_download", 0.0 },
        { "recon_bcq", 0.0 },
        { "recon_kernel", 0.0 },
        { "recon_limiter_upload", 0.0 },
        { "recon_limiter_kernel", 0.0 },
        { "bnd_refresh", 0.0 },
        { "bnd_static_upload", 0.0 },
        { "bnd_kernel", 0.0 },
        { "bnd_ghost_copy", 0.0 },
        { "bnd_geometry", 0.0 },
        { "bnd_prep", 0.0 },
        { "bnd_backend", 0.0 },
        { "bnd_resize", 0.0 },
        // Macro-step segments. Record() matches by name against this table and
        // silently drops anything that is not listed here, so a new timer has
        // to be registered in both places or it never shows up in the TSV.
        { "macro_time_step", 0.0 },
        { "macro_update_residuals", 0.0 },
        { "macro_calc_lhs", 0.0 },
        { "macro_update_flowfield", 0.0 },
        { "macro_calc_boundary", 0.0 },
        { "macro_interface", 0.0 },
        { "g_prep", 0.0 },
        { "g_contract", 0.0 },
        { "g_geom", 0.0 },
        { "v_upload", 0.0 },
        { "v_kernel", 0.0 },
        { "v_download", 0.0 },
        { "r_prep", 0.0 },
        { "fused_scale", 0.0 },
    };
    return totals;
}

// Record() only accumulates into the fixed table above.  A timer whose name is
// missing from that table used to disappear without a trace, which made a whole
// instrumentation round look like "that segment costs nothing".  Count the
// misses and report them on Flush so the mistake is visible next time.
struct DroppedCategory
{
    unsigned long long count = 0;
    std::string first;
};

DroppedCategory & DroppedCategories()
{
    static DroppedCategory dropped;
    return dropped;
}

std::string RankAwarePath( const char * configured )
{
    std::string path( configured == nullptr ? "" : configured );
    const std::string token( "%r" );
    const std::size_t tokenPos = path.find( token );
    if ( tokenPos == std::string::npos ) return path;

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
    path.replace( tokenPos, token.size(), rank );
    return path;
}

}

bool StageProfiler::Enabled()
{
    // Called from every ScopedStageTimer construction, i.e. hundreds of times
    // per RK step.  The switch is a process-lifetime decision, so resolve it
    // once instead of scanning the environment each time.
    static const bool enabled = []()
    {
        const char * path = std::getenv( "ONEFLOW_STAGE_BREAKDOWN_FILE" );
        return path != nullptr && path[ 0 ] != '\0';
    }();
    return enabled;
}

void StageProfiler::Record( const char * category, double elapsedMs ) noexcept
{
    if ( category == nullptr || category[ 0 ] == '\0' || elapsedMs < 0.0 )
        return;

    std::lock_guard< std::mutex > lock( OutputMutex() );
    for ( CategoryTotal & total : Totals() )
    {
        // Unused trailing slots keep a null name; never compare against them.
        if ( total.name == nullptr ) continue;
        if ( std::strcmp( total.name, category ) == 0 )
        {
            total.elapsedMs += elapsedMs;
            return;
        }
    }
    DroppedCategory & dropped = DroppedCategories();
    ++dropped.count;
    try
    {
        if ( dropped.first.empty() ) dropped.first = category;
    }
    catch ( ... )
    {
        // Profiling must never change solver behaviour.
    }
}

void StageProfiler::RecordTask( int taskId, const char * taskName,
                                double elapsedMs ) noexcept
{
    if ( elapsedMs < 0.0 ) return;
    try
    {
        std::lock_guard< std::mutex > lock( OutputMutex() );
        auto & totals = TaskTotals();
        auto found = totals.find( taskId );
        if ( found == totals.end() )
        {
            TaskTotal total;
            total.name = taskName == nullptr
                ? std::string( "<unnamed>" ) : std::string( taskName );
            found = totals.emplace( taskId, std::move( total ) ).first;
        }
        found->second.count += 1;
        found->second.elapsedMs += elapsedMs;
    }
    catch ( ... )
    {
        // Profiling must never change solver behaviour.
    }
}

void StageProfiler::Flush()
{
    const char * configured = std::getenv( "ONEFLOW_STAGE_BREAKDOWN_FILE" );
    if ( configured == nullptr || configured[ 0 ] == '\0' ) return;

    const std::string path = RankAwarePath( configured );
    if ( path.empty() ) return;

    std::lock_guard< std::mutex > lock( OutputMutex() );
    std::ofstream output( path, std::ios::out | std::ios::trunc );
    if ( ! output )
    {
        throw std::runtime_error(
            "cannot open ONEFLOW_STAGE_BREAKDOWN_FILE: " + path );
    }
    output << "category\telapsed_ms\n";
    for ( const CategoryTotal & total : Totals() )
    {
        if ( total.name == nullptr ) continue;
        output << total.name << "\t" << total.elapsedMs << "\n";
    }
    if ( DroppedCategories().count != 0 )
    {
        std::cerr << "StageProfiler: dropped " << DroppedCategories().count
                  << " timer records for unregistered category \""
                  << DroppedCategories().first
                  << "\"; register it in Totals()" << std::endl;
    }
    if ( ! output )
    {
        throw std::runtime_error(
            "failed while writing ONEFLOW_STAGE_BREAKDOWN_FILE: " + path );
    }

    const std::string taskPath = path + ".tasks.tsv";
    std::vector< std::pair< int, const TaskTotal * > > ordered;
    ordered.reserve( TaskTotals().size() );
    for ( const auto & entry : TaskTotals() )
    {
        ordered.emplace_back( entry.first, & entry.second );
    }
    std::sort( ordered.begin(), ordered.end(),
        []( const std::pair< int, const TaskTotal * > & left,
            const std::pair< int, const TaskTotal * > & right )
        {
            return left.second->elapsedMs > right.second->elapsedMs;
        } );

    std::ofstream taskOutput( taskPath, std::ios::out | std::ios::trunc );
    if ( ! taskOutput )
    {
        throw std::runtime_error(
            "cannot open task breakdown file: " + taskPath );
    }
    taskOutput << "task_id\ttask\tcount\ttotal_ms\tmean_ms\n";
    for ( const auto & entry : ordered )
    {
        const TaskTotal & total = * entry.second;
        taskOutput << entry.first << '\t' << total.name << '\t'
                   << total.count << '\t' << total.elapsedMs << '\t'
                   << ( total.count == 0 ? 0.0
                                         : total.elapsedMs / total.count )
                   << '\n';
    }
    if ( ! taskOutput )
    {
        throw std::runtime_error(
            "failed while writing task breakdown file: " + taskPath );
    }

    if ( HostAccountingEnabled() )
    {
        if ( DroppedHostRecords() != 0 )
        {
            throw std::runtime_error( "host accounting lost timer records" );
        }
        std::ofstream hostOutput( path + ".host.tsv" );
        hostOutput << std::setprecision( 17 );
        hostOutput << "scope\tcount\tinclusive_ms\texclusive_ms\n";
        for ( const auto & entry : HostTotals() )
        {
            const HostTotal & total = entry.second;
            hostOutput << entry.first << '\t' << total.count << '\t'
                       << total.inclusiveMs << '\t' << total.exclusiveMs << '\n';
        }
        if ( ! hostOutput )
        {
            throw std::runtime_error( "cannot write host accounting file: " + path );
        }
    }
}

void StageProfiler::RecordHost( const char * scope, double inclusiveMs,
                                double exclusiveMs ) noexcept
{
    std::lock_guard< std::mutex > lock( OutputMutex() );
    try
    {
        HostTotal & total = HostTotals()[ scope ];
        ++total.count;
        total.inclusiveMs += inclusiveMs;
        total.exclusiveMs += exclusiveMs;
    }
    catch ( ... )
    {
        ++DroppedHostRecords();
    }
}

thread_local ScopedHostTimer * ScopedHostTimer::current_ = nullptr;

ScopedHostTimer::ScopedHostTimer( const char * scope, bool root )
    : scope_( scope )
{
    if ( ! HostAccountingEnabled() || scope == nullptr ) return;
    // Initialization and final output are outside the steady solver root.
    if ( current_ == nullptr && ! root ) return;
    parent_ = current_;
    active_ = true;
    start_ = std::chrono::steady_clock::now();
    current_ = this;
}

ScopedHostTimer::~ScopedHostTimer() noexcept
{
    if ( ! active_ ) return;
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration< double, std::milli > elapsed = end - start_;
    current_ = parent_;
    if ( parent_ != nullptr ) parent_->childrenMs_ += elapsed.count();
    StageProfiler::RecordHost( scope_, elapsed.count(),
                              std::max( 0.0, elapsed.count() - childrenMs_ ) );
}

const char * StageProfiler::CategoryForTask( const std::string & taskName )
{
    if ( taskName == "LOAD_Q" ) return "load_q";
    if ( taskName == "CALC_TIME_STEP" ) return "calc_time_step";
    if ( taskName == "LOAD_RESIDUALS" ) return "load_residuals";
    if ( taskName == "UPDATE_RESIDUALS" ) return "update_residuals";
    if ( taskName == "CALC_LHS" ) return "calc_lhs";
    if ( taskName == "UPDATE_FLOWFIELD" ) return "update_flowfield";
    if ( taskName == "CALC_BOUNDARY" ) return "calc_boundary";
    if ( taskName.find( "INIT" ) != std::string::npos
         || taskName.find( "READ_RESTART" ) != std::string::npos )
    {
        return "initialization";
    }
    if ( taskName.find( "INTERFACE" ) != std::string::npos
         || taskName.find( "EXCHANGE" ) != std::string::npos
         || taskName.find( "UPLOAD_INTERFACE" ) != std::string::npos
         || taskName.find( "DOWNLOAD_INTERFACE" ) != std::string::npos )
    {
        return "mpi_interface";
    }
    if ( taskName.find( "DUMP" ) != std::string::npos
         || taskName.find( "WRITE" ) != std::string::npos
         || taskName.find( "VISUALIZATION" ) != std::string::npos )
    {
        return "output";
    }
    return "";
}

ScopedStageTimer::ScopedStageTimer( const char * category )
    : category_( category ), enabled_( StageProfiler::Enabled() )
{
    if ( enabled_ )
    {
        start_ = std::chrono::steady_clock::now();
    }
}

ScopedStageTimer::~ScopedStageTimer() noexcept
{
    if ( ! enabled_ ) return;
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration< double, std::milli > elapsed = end - start_;
    StageProfiler::Record( category_, elapsed.count() );
}

EndNameSpace
