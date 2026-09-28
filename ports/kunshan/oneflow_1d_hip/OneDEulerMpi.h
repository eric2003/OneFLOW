#pragma once

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

namespace oneflow_1d
{

struct EulerMpiProblem
{
    int globalNx = 0;
    int localNx = 0;
    int globalOffset = 0;
    int steps = 0;
    double gamma = 1.4;
    double dt = 0.0;
    double dx = 0.0;
    MPI_Comm communicator = MPI_COMM_NULL;
};

struct EulerMpiTiming
{
    double createMilliseconds = 0.0;
    double uploadMilliseconds = 0.0;
    double advanceMilliseconds = 0.0;
    double downloadMilliseconds = 0.0;
    double mpiMilliseconds = 0.0;
    double computeMilliseconds = 0.0;
    double kernelMilliseconds = 0.0;
    int kernelLaunches = 0;
    int haloExchanges = 0;
    int deviceSynchronizations = 0;
};

struct EulerMpiResult
{
    std::vector< double > state;
    EulerMpiTiming timing;
    std::string deviceName;
    std::string architecture;
    int localRank = -1;
    int localRanks = 0;
    int deviceIndex = -1;
    int visibleDeviceCount = 0;
};

EulerMpiResult RunCpuEulerMpi(
    const EulerMpiProblem & problem,
    const std::vector< double > & initialState );

#ifdef ONEFLOW_1D_MPI_USE_HIP
EulerMpiResult RunHipEulerMpi(
    const EulerMpiProblem & problem,
    const std::vector< double > & initialState );
#endif

std::uint64_t EulerMpiStateHash(
    const std::vector< double > & localState, int globalNx, int localNx,
    int globalOffset );

} // namespace oneflow_1d
