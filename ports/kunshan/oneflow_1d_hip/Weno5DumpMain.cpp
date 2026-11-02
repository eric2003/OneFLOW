#include "OneDWeno5.h"
#include "Weno5DumpFormat.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
constexpr double Gamma = 1.4;

std::vector<double> SmoothState( int nx )
{
    std::vector<double> state( oneflow_1d::EulerComponents * nx );
    for ( int i = 0; i < nx; ++ i )
    {
        const double x = ( i + 0.5 ) / nx;
        const double rho = 1.0 + 0.1 * std::sin( 2.0 * M_PI * x );
        const double velocity = 0.2 + 0.05 * std::cos( 2.0 * M_PI * x );
        const double pressure = 1.0 + 0.08 * std::sin( 4.0 * M_PI * x );
        state[i] = rho;
        state[nx + i] = rho * velocity;
        state[2 * nx + i] = pressure / ( Gamma - 1.0 )
            + 0.5 * rho * velocity * velocity;
    }
    return state;
}

std::vector<double> SodState( int nx )
{
    std::vector<double> state( oneflow_1d::EulerComponents * nx );
    for ( int i = 0; i < nx; ++ i )
    {
        const double x = ( i + 0.5 ) / nx;
        const double rho = x <= 0.5 ? 1.0 : 0.125;
        const double pressure = x <= 0.5 ? 1.0 : 0.1;
        state[i] = rho;
        state[nx + i] = 0.0;
        state[2 * nx + i] = pressure / ( Gamma - 1.0 );
    }
    return state;
}

oneflow_1d::Weno5CaseHeader Header( const char * name, int nx, int steps,
    oneflow_1d::EulerBoundary boundary, double dt )
{
    oneflow_1d::Weno5CaseHeader header;
    std::strncpy( header.name.data(), name, header.name.size() - 1 );
    header.nx = nx;
    header.steps = steps;
    header.boundary = static_cast<std::int32_t>( boundary );
    header.gamma = Gamma;
    header.dt = dt;
    header.dx = 1.0 / nx;
    return header;
}

void WriteCase( std::ostream & output, const oneflow_1d::Weno5CaseHeader & header,
                std::vector<double> state )
{
    oneflow_1d::WritePod( output, header );
    oneflow_1d::WriteValues( output, state );
    for ( int step = 0; step < header.steps; ++ step )
    {
        oneflow_1d::Weno5Trace trace;
        oneflow_1d::OneDCpuLaxWeno5Step( state.data(), header.nx, header.gamma,
            header.dt, header.dx,
            static_cast<oneflow_1d::EulerBoundary>( header.boundary ), trace );
        oneflow_1d::WriteTrace( output, trace );
        const int cellValues = oneflow_1d::EulerComponents * header.nx;
        std::copy( trace.state.end() - cellValues, trace.state.end(), state.begin() );
    }
}
} // namespace

int main( int argc, char ** argv )
{
    try
    {
        const std::string path = argc > 1 ? argv[1] : "weno5_cpu_trace.bin";
        std::ofstream output( path, std::ios::binary | std::ios::trunc );
        if ( ! output ) throw std::runtime_error( "cannot open dump output: " + path );
        output.write( oneflow_1d::Weno5DumpMagic.data(), oneflow_1d::Weno5DumpMagic.size() );
        oneflow_1d::WritePod( output, oneflow_1d::Weno5DumpVersion );
        constexpr std::uint32_t caseCount = 2;
        oneflow_1d::WritePod( output, caseCount );
        constexpr int nx = 257;
        constexpr int steps = 10;
        WriteCase( output, Header( "smooth-periodic", nx, steps,
            oneflow_1d::EulerBoundary::Periodic, 0.05 / nx ), SmoothState( nx ) );
        WriteCase( output, Header( "sod-transmissive", nx, steps,
            oneflow_1d::EulerBoundary::Transmissive, 1.0e-4 ), SodState( nx ) );
        output.close();
        if ( ! output ) throw std::runtime_error( "failed to close dump output" );
        std::cout << "WENO5 CPU reference dump: " << path << "\n";
        std::cout << "cases=2 nx=" << nx << " steps=" << steps
                  << " components=3 format_version=" << oneflow_1d::Weno5DumpVersion << "\n";
        return 0;
    }
    catch ( const std::exception & error )
    {
        std::cerr << "WENO5 CPU dump failed: " << error.what() << "\n";
        return 1;
    }
}
