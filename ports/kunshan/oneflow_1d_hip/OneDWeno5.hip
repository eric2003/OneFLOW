#include "OneDWeno5.h"

#include <hip/hip_runtime.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace oneflow_1d
{
namespace
{
void Check( hipError_t error, const char * operation )
{
    if ( error != hipSuccess )
        throw std::runtime_error( std::string( operation ) + ": " + hipGetErrorString( error ) );
}

__device__ int PI( int c, int i, int nx )
{
    return c * ( nx + 2 * WenoGhostCells ) + i + WenoGhostCells;
}
__device__ int CI( int c, int i, int nx ) { return c * nx + i; }
__device__ int FI( int c, int i, int nx ) { return c * ( nx + 1 ) + i; }
__device__ double Sq( double x ) { return x * x; }

__device__ double Wl( double v1, double v2, double v3, double v4, double v5 )
{
    constexpr double eps = 1.0e-6;
    const double s1 = 13.0 / 12.0 * Sq( v1 - 2.0 * v2 + v3 )
        + 0.25 * Sq( v1 - 4.0 * v2 + 3.0 * v3 );
    const double s2 = 13.0 / 12.0 * Sq( v2 - 2.0 * v3 + v4 )
        + 0.25 * Sq( v2 - v4 );
    const double s3 = 13.0 / 12.0 * Sq( v3 - 2.0 * v4 + v5 )
        + 0.25 * Sq( 3.0 * v3 - 4.0 * v4 + v5 );
    const double a1 = 0.1 / Sq( eps + s1 );
    const double a2 = 0.6 / Sq( eps + s2 );
    const double a3 = 0.3 / Sq( eps + s3 );
    const double sum = a1 + a2 + a3;
    const double q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3;
    const double q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0;
    const double q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0;
    return a1 / sum * q1 + a2 / sum * q2 + a3 / sum * q3;
}

__device__ double Wr( double v1, double v2, double v3, double v4, double v5 )
{
    constexpr double eps = 1.0e-6;
    const double s1 = 13.0 / 12.0 * Sq( v1 - 2.0 * v2 + v3 )
        + 0.25 * Sq( v1 - 4.0 * v2 + 3.0 * v3 );
    const double s2 = 13.0 / 12.0 * Sq( v2 - 2.0 * v3 + v4 )
        + 0.25 * Sq( v2 - v4 );
    const double s3 = 13.0 / 12.0 * Sq( v3 - 2.0 * v4 + v5 )
        + 0.25 * Sq( 3.0 * v3 - 4.0 * v4 + v5 );
    const double a1 = 0.3 / Sq( eps + s1 );
    const double a2 = 0.6 / Sq( eps + s2 );
    const double a3 = 0.1 / Sq( eps + s3 );
    const double sum = a1 + a2 + a3;
    const double q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + v3 / 3.0;
    const double q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0;
    const double q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0;
    return a1 / sum * q1 + a2 / sum * q2 + a3 / sum * q3;
}

__global__ void FillPadded( const double * state, int nx, int boundary, double * padded )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = nx + 2 * WenoGhostCells;
    if ( index >= EulerComponents * stride ) return;
    const int component = index / stride;
    const int cell = index - component * stride - WenoGhostCells;
    if ( cell >= 0 && cell < nx )
        padded[index] = state[CI( component, cell, nx )];
    else if ( boundary == static_cast<int>( EulerBoundary::Periodic ) )
        padded[index] = state[CI( component, cell < 0 ? nx + cell : cell - nx, nx )];
    else
        padded[index] = state[CI( component, cell < 0 ? 0 : nx - 1, nx )];
}

__global__ void SplitFlux( const double * padded, int nx, double gamma,
                           double * positive, double * negative )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = nx + 2 * WenoGhostCells;
    if ( index >= EulerComponents * stride ) return;
    const int component = index / stride;
    const int local = index - component * stride;
    const double rho = padded[local];
    const double momentum = padded[stride + local];
    const double energy = padded[2 * stride + local];
    const double velocity = momentum / rho;
    const double specificEnergy = energy / rho;
    const double pressure = ( gamma - 1.0 )
        * ( rho * specificEnergy - 0.5 * rho * velocity * velocity );
    const double enthalpy = specificEnergy + pressure / rho;
    const double weightedVelocity = sqrt( fabs( rho ) ) * velocity;
    const double weightedEnthalpy = sqrt( fabs( rho ) ) * enthalpy;
    const double soundSpeed = sqrt( fabs( ( gamma - 1.0 )
        * ( weightedEnthalpy - 0.5 * weightedVelocity * weightedVelocity ) ) );
    const double speed = fabs( soundSpeed + weightedVelocity );
    const double physical = component == 0 ? momentum
        : ( component == 1 ? momentum * velocity + pressure
                           : momentum * energy / rho + pressure * momentum / rho );
    positive[index] = 0.5 * ( physical + speed * padded[index] );
    negative[index] = 0.5 * ( physical - speed * padded[index] );
}

__global__ void Reconstruct( const double * positive, const double * negative,
                             int nx, double * left, double * right, double * flux )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = nx + 1;
    if ( index >= EulerComponents * stride ) return;
    const int component = index / stride;
    const int face = index - component * stride;
    const int center = face - 1;
    left[index] = Wl(
        positive[PI( component, center - 2, nx )],
        positive[PI( component, center - 1, nx )],
        positive[PI( component, center, nx )],
        positive[PI( component, center + 1, nx )],
        positive[PI( component, center + 2, nx )] );
    right[index] = Wr(
        negative[PI( component, center - 1, nx )],
        negative[PI( component, center, nx )],
        negative[PI( component, center + 1, nx )],
        negative[PI( component, center + 2, nx )],
        negative[PI( component, center + 3, nx )] );
    flux[index] = left[index] + right[index];
}

__global__ void AssembleResidual( const double * positive, const double * negative,
                                  int nx, double dx, double * residual )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index >= EulerComponents * nx ) return;
    const int component = index / nx;
    const int cell = index - component * nx;
    residual[index] = -(
        ( positive[FI( component, cell + 1, nx )]
            - positive[FI( component, cell, nx )] ) / dx
        + ( negative[FI( component, cell + 1, nx )]
            - negative[FI( component, cell, nx )] ) / dx );
}

__global__ void Update( const double * base, const double * current,
                        const double * residual, int count, double dt,
                        double baseWeight, double updateWeight, double * next )
{
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if ( index < count ) next[index] = baseWeight * base[index]
        + updateWeight * ( current[index] + dt * residual[index] );
}

} // namespace

void OneDHipLaxWeno5Step( const double * input, int nx, double gamma, double dt,
                          double dx, EulerBoundary boundary, Weno5Trace & trace )
{
    if ( input == nullptr || nx < WenoGhostCells || gamma <= 1.0
        || dt <= 0.0 || dx <= 0.0 )
        throw std::invalid_argument( "invalid HIP Lax-WENO5 input" );
    ResizeWeno5Trace( trace, nx );
    const int paddedValues = EulerComponents * ( nx + 2 * WenoGhostCells );
    const int faceValues = EulerComponents * ( nx + 1 );
    const int cellValues = EulerComponents * nx;
    constexpr int threads = 256;
    const int paddedBlocks = ( paddedValues + threads - 1 ) / threads;
    const int faceBlocks = ( faceValues + threads - 1 ) / threads;
    const int cellBlocks = ( cellValues + threads - 1 ) / threads;
    const size_t paddedBytes = paddedValues * sizeof( double );
    const size_t faceBytes = faceValues * sizeof( double );
    const size_t cellBytes = cellValues * sizeof( double );
    double * base = nullptr, * current = nullptr, * next = nullptr;
    double * padded = nullptr, * positive = nullptr, * negative = nullptr;
    double * left = nullptr, * right = nullptr, * flux = nullptr, * residual = nullptr;
    auto cleanup = [&]() {
        hipFree( residual ); hipFree( flux ); hipFree( right ); hipFree( left );
        hipFree( negative ); hipFree( positive ); hipFree( padded );
        hipFree( next ); hipFree( current ); hipFree( base );
    };
    try
    {
        Check( hipMalloc( &base, cellBytes ), "malloc base" );
        Check( hipMalloc( &current, cellBytes ), "malloc current" );
        Check( hipMalloc( &next, cellBytes ), "malloc next" );
        Check( hipMalloc( &padded, paddedBytes ), "malloc padded" );
        Check( hipMalloc( &positive, paddedBytes ), "malloc split positive" );
        Check( hipMalloc( &negative, paddedBytes ), "malloc split negative" );
        Check( hipMalloc( &left, faceBytes ), "malloc reconstructed positive" );
        Check( hipMalloc( &right, faceBytes ), "malloc reconstructed negative" );
        Check( hipMalloc( &flux, faceBytes ), "malloc numerical flux" );
        Check( hipMalloc( &residual, cellBytes ), "malloc residual" );
        Check( hipMemcpy( base, input, cellBytes, hipMemcpyHostToDevice ), "copy base" );
        Check( hipMemcpy( current, input, cellBytes, hipMemcpyHostToDevice ), "copy current" );
        std::copy( input, input + cellValues, trace.state.begin() );
        for ( int stage = 0; stage < EulerRkStages; ++ stage )
        {
            hipLaunchKernelGGL( FillPadded, dim3( paddedBlocks ), dim3( threads ),
                0, 0, current, nx, static_cast<int>( boundary ), padded );
            Check( hipGetLastError(), "fill padded launch" );
            hipLaunchKernelGGL( SplitFlux, dim3( paddedBlocks ), dim3( threads ),
                0, 0, padded, nx, gamma, positive, negative );
            Check( hipGetLastError(), "split flux launch" );
            hipLaunchKernelGGL( Reconstruct, dim3( faceBlocks ), dim3( threads ),
                0, 0, positive, negative, nx, left, right, flux );
            Check( hipGetLastError(), "reconstruct launch" );
            hipLaunchKernelGGL( AssembleResidual, dim3( cellBlocks ), dim3( threads ),
                0, 0, left, right, nx, dx, residual );
            Check( hipGetLastError(), "residual launch" );
            const double baseWeight = stage == 1 ? 0.75
                : ( stage == 2 ? 1.0 / 3.0 : 0.0 );
            const double updateWeight = stage == 0 ? 1.0
                : ( stage == 1 ? 0.25 : 2.0 / 3.0 );
            hipLaunchKernelGGL( Update, dim3( cellBlocks ), dim3( threads ),
                0, 0, base, current, residual, cellValues, dt,
                baseWeight, updateWeight, next );
            Check( hipGetLastError(), "update launch" );
            Check( hipDeviceSynchronize(), "stage synchronize" );
            const int paddedOffset = stage * paddedValues;
            const int faceOffset = stage * faceValues;
            const int cellOffset = stage * cellValues;
            Check( hipMemcpy( trace.splitPositive.data() + paddedOffset, positive,
                paddedBytes, hipMemcpyDeviceToHost ), "copy split positive" );
            Check( hipMemcpy( trace.splitNegative.data() + paddedOffset, negative,
                paddedBytes, hipMemcpyDeviceToHost ), "copy split negative" );
            Check( hipMemcpy( trace.reconstructedPositive.data() + faceOffset, left,
                faceBytes, hipMemcpyDeviceToHost ), "copy reconstructed positive" );
            Check( hipMemcpy( trace.reconstructedNegative.data() + faceOffset, right,
                faceBytes, hipMemcpyDeviceToHost ), "copy reconstructed negative" );
            Check( hipMemcpy( trace.numericalFlux.data() + faceOffset, flux,
                faceBytes, hipMemcpyDeviceToHost ), "copy numerical flux" );
            Check( hipMemcpy( trace.residual.data() + cellOffset, residual,
                cellBytes, hipMemcpyDeviceToHost ), "copy residual" );
            Check( hipMemcpy( trace.state.data() + ( stage + 1 ) * cellValues, next,
                cellBytes, hipMemcpyDeviceToHost ), "copy state" );
            std::swap( current, next );
        }
    }
    catch ( ... ) { cleanup(); throw; }
    cleanup();
}

} // namespace oneflow_1d
