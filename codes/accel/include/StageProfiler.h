/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
-------------------------------------------------------------------------------
    Opt-in stage timing for accelerator/dataflow investigations.

    The profiler is deliberately inert unless ONEFLOW_STAGE_BREAKDOWN_FILE is
    set.  Each MPI process should normally use a path containing "%r" so that
    independent task records cannot interleave.
\*---------------------------------------------------------------------------*/

#pragma once

#include "NamespaceMacros.h"

#include <chrono>
#include <string>

BeginNameSpace( ONEFLOW )

class StageProfiler
{
public:
    static bool Enabled();

    // Accumulate one completed event in process memory.  Flush writes one
    // total per category so profiling does not add file I/O to every stage.
    static void Record( const char * category, double elapsedMs ) noexcept;

    // Write the accumulated totals to the configured TSV.  This is called at
    // the end of the solver pipeline so output failures propagate as workload
    // failures instead of being hidden by process teardown.
    static void Flush();

    // Map task names to the coarse categories requested by the benchmark
    // contract.  An empty string means that the task is not a coarse-stage
    // boundary and should not be timed here.
    static const char * CategoryForTask( const std::string & taskName );

    // Accumulate the end-to-end dispatch time of one task (plan + execute).
    // Unlike Record(), this is keyed by operation id so the task breakdown can
    // show which dispatch owns the time that the coarse categories do not
    // name.  Flush() writes these totals to "<breakdown>.tasks.tsv".
    static void RecordTask( int taskId, const char * taskName,
                            double elapsedMs ) noexcept;
};

class ScopedStageTimer
{
public:
    explicit ScopedStageTimer( const char * category );
    ~ScopedStageTimer() noexcept;

    ScopedStageTimer( const ScopedStageTimer & ) = delete;
    ScopedStageTimer & operator=( const ScopedStageTimer & ) = delete;

private:
    const char * category_ = nullptr;
    bool enabled_ = false;
    std::chrono::steady_clock::time_point start_;
};

EndNameSpace
