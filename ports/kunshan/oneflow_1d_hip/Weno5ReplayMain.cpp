#include "OneDWeno5.h"
#include "Weno5DumpFormat.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
struct Metrics
{
    double abs = 0.0;
    double scaledRelative = 0.0;
    std::uint64_t ulp = 0;
    bool ok = true;
};

std::uint64_t OrderedBits( double value )
{
    std::uint64_t bits = 0;
    std::memcpy( &bits, &value, sizeof( bits ) );
    return ( bits >> 63 ) != 0 ? ~bits : bits | ( std::uint64_t( 1 ) << 63 );
}

void Compare( Metrics & metrics, double reference, double candidate )
{
    if ( !std::isfinite( reference ) || !std::isfinite( candidate ) )
    {
        metrics.ok = false;
        return;
    }
    const double error = std::abs( candidate - reference );
    const double scale = std::max( 1.0, std::abs( reference ) );
    metrics.abs = std::max( metrics.abs, error );
    metrics.scaledRelative = std::max( metrics.scaledRelative, error / scale );
    const std::uint64_t a = OrderedBits( reference );
    const std::uint64_t b = OrderedBits( candidate );
    metrics.ulp = std::max( metrics.ulp, a > b ? a - b : b - a );
    if ( error > 1.0e-15 + 1.0e-15 * scale ) metrics.ok = false;
}

bool SmoothCell( const std::vector<double> & state, int stage, int cell, int nx, double gamma )
{
    if ( cell <= 0 || cell >= nx - 1 ) return false;
    const int offset = stage * oneflow_1d::EulerComponents * nx;
    for ( int component : {0, 2} )
    {
        const auto value = [&]( int i ) {
            const double rho = state[offset + i];
            const double momentum = state[offset + nx + i];
            const double energy = state[offset + 2 * nx + i];
            if ( component == 0 ) return rho;
            return ( gamma - 1.0 ) * ( energy - 0.5 * momentum * momentum / rho );
        };
        const double center = value( cell );
        const double scale = std::max( 1.0, std::abs( center ) );
        if ( std::abs( center - value( cell - 1 ) ) > 0.02 * scale
            || std::abs( center - value( cell + 1 ) ) > 0.02 * scale ) return false;
    }
    return true;
}

bool Physical( const std::vector<double> & state, int nx, double gamma )
{
    const int cellValues = oneflow_1d::EulerComponents * nx;
    for ( int stage = 0; stage <= oneflow_1d::EulerRkStages; ++ stage )
    for ( int i = 0; i < nx; ++ i )
    {
        const int offset = stage * cellValues;
        const double rho = state[offset + i];
        const double momentum = state[offset + nx + i];
        const double energy = state[offset + 2 * nx + i];
        const double pressure = ( gamma - 1.0 )
            * ( energy - 0.5 * momentum * momentum / rho );
        if ( !std::isfinite( rho ) || !std::isfinite( pressure )
            || rho <= 0.0 || pressure <= 0.0 ) return false;
    }
    return true;
}

struct CaseMetrics
{
    std::array<Metrics, 3> splitPositive;
    std::array<Metrics, 3> splitNegative;
    std::array<Metrics, 3> reconstructedPositive;
    std::array<Metrics, 3> reconstructedNegative;
    std::array<Metrics, 3> numericalFlux;
    std::array<Metrics, 3> residual;
    std::array<Metrics, 3> state;
    std::array<Metrics, 3> primitive;
    int excludedResidual = 0;
    bool physical = true;
};

void CompareResidual( std::array<Metrics, 3> & metrics,
    const std::vector<double> & reference, const std::vector<double> & candidate,
    int offset, bool useSmoothMask, const std::vector<double> & referenceState,
    int stage, int nx, double gamma, int & excluded )
{
    for ( int component = 0; component < oneflow_1d::EulerComponents; ++ component )
    {
        for ( int cell = 0; cell < nx; ++ cell )
        {
            if ( useSmoothMask && !SmoothCell( referenceState, stage, cell, nx, gamma ) )
            {
                if ( component == 0 ) ++ excluded;
                continue;
            }
            const int index = offset + component * nx + cell;
            Compare( metrics[component], reference[index], candidate[index] );
        }
    }
}

void Report( const char * label, const std::array<Metrics, 3> & metrics )
{
    for ( int component = 0; component < 3; ++ component )
    {
        const Metrics & m = metrics[component];
        std::cout << "  " << label << " U" << component
                  << " abs=" << m.abs << " rel=" << m.scaledRelative
                  << " ulp=" << m.ulp << " " << ( m.ok ? "PASS" : "FAIL" ) << "\n";
    }
}

bool RunCase( std::istream & input, const oneflow_1d::Weno5CaseHeader & header )
{
    const int nx = header.nx;
    const int pv = oneflow_1d::EulerComponents * ( nx + 2 * oneflow_1d::WenoGhostCells );
    const int fv = oneflow_1d::EulerComponents * ( nx + 1 );
    const int cv = oneflow_1d::EulerComponents * nx;
    std::vector<double> state( cv );
    oneflow_1d::ReadValues( input, state );
    oneflow_1d::Weno5Trace reference;
    bool pass = true;
    CaseMetrics metrics;
    for ( int step = 0; step < header.steps; ++ step )
    {
        oneflow_1d::ReadTrace( input, nx, reference );
        oneflow_1d::Weno5Trace candidate;
        oneflow_1d::OneDHipLaxWeno5Step( state.data(), nx, header.gamma,
            header.dt, header.dx, static_cast<oneflow_1d::EulerBoundary>( header.boundary ), candidate );
        for ( int stage = 0; stage < oneflow_1d::EulerRkStages; ++ stage )
        {
            const int po = stage * pv;
            const int fo = stage * fv;
            const int co = stage * cv;
            for ( int component = 0; component < 3; ++ component )
            {
                for ( int i = 0; i < nx + 2 * oneflow_1d::WenoGhostCells; ++ i )
                {
                    Compare( metrics.splitPositive[component], reference.splitPositive[po + component * ( nx + 2 * oneflow_1d::WenoGhostCells ) + i], candidate.splitPositive[po + component * ( nx + 2 * oneflow_1d::WenoGhostCells ) + i] );
                    Compare( metrics.splitNegative[component], reference.splitNegative[po + component * ( nx + 2 * oneflow_1d::WenoGhostCells ) + i], candidate.splitNegative[po + component * ( nx + 2 * oneflow_1d::WenoGhostCells ) + i] );
                }
                for ( int i = 0; i < nx + 1; ++ i )
                {
                    Compare( metrics.reconstructedPositive[component], reference.reconstructedPositive[fo + component * ( nx + 1 ) + i], candidate.reconstructedPositive[fo + component * ( nx + 1 ) + i] );
                    Compare( metrics.reconstructedNegative[component], reference.reconstructedNegative[fo + component * ( nx + 1 ) + i], candidate.reconstructedNegative[fo + component * ( nx + 1 ) + i] );
                    Compare( metrics.numericalFlux[component], reference.numericalFlux[fo + component * ( nx + 1 ) + i], candidate.numericalFlux[fo + component * ( nx + 1 ) + i] );
                }
            }
            CompareResidual( metrics.residual, reference.residual, candidate.residual, co,
                header.boundary == static_cast<int>( oneflow_1d::EulerBoundary::Transmissive ),
                reference.state, stage, nx, header.gamma, metrics.excludedResidual );
            for ( int component = 0; component < 3; ++ component )
                for ( int i = 0; i < nx; ++ i )
                    Compare( metrics.state[component], reference.state[( stage + 1 ) * cv + component * nx + i], candidate.state[( stage + 1 ) * cv + component * nx + i] );
            for ( int i = 0; i < nx; ++ i )
            {
                const auto primitive = [&]( const std::vector<double> & values, int component ) {
                    const double rho = values[( stage + 1 ) * cv + i];
                    const double momentum = values[( stage + 1 ) * cv + nx + i];
                    const double energy = values[( stage + 1 ) * cv + 2 * nx + i];
                    if ( component == 0 ) return rho;
                    if ( component == 1 ) return momentum / rho;
                    return ( header.gamma - 1.0 ) * ( energy - 0.5 * momentum * momentum / rho );
                };
                Compare( metrics.primitive[0], primitive( reference.state, 0 ), primitive( candidate.state, 0 ) );
                Compare( metrics.primitive[1], primitive( reference.state, 1 ), primitive( candidate.state, 1 ) );
                Compare( metrics.primitive[2], primitive( reference.state, 2 ), primitive( candidate.state, 2 ) );
            }
        }
        metrics.physical = metrics.physical
            && Physical( reference.state, nx, header.gamma )
            && Physical( candidate.state, nx, header.gamma );
        std::copy( candidate.state.end() - cv, candidate.state.end(), state.begin() );
    }
    auto all = [&]( const auto & values ) { for ( const auto & value : values ) pass = pass && value.ok; };
    all( metrics.splitPositive ); all( metrics.splitNegative ); all( metrics.reconstructedPositive ); all( metrics.reconstructedNegative );
    all( metrics.numericalFlux ); all( metrics.residual ); all( metrics.state ); all( metrics.primitive );
    pass = pass && metrics.physical;
    std::cout << header.name.data() << " steps=" << header.steps
              << " residual_smooth_excluded=" << metrics.excludedResidual << "\n";
    Report( "split+", metrics.splitPositive ); Report( "split-", metrics.splitNegative );
    Report( "recon+", metrics.reconstructedPositive ); Report( "recon-", metrics.reconstructedNegative );
    Report( "flux", metrics.numericalFlux ); Report( "residual", metrics.residual );
    Report( "state", metrics.state ); Report( "primitive", metrics.primitive );
    std::cout << "  physical_state=" << ( metrics.physical ? "PASS" : "FAIL" )
              << " result=" << ( pass ? "PASS" : "FAIL" ) << "\n";
    return pass;
}
}

int main( int argc, char ** argv )
{
    std::cout << std::setprecision( 17 ) << std::scientific;
    try
    {
        const std::string path = argc > 1 ? argv[1] : "weno5_cpu_trace.bin";
        std::ifstream input( path, std::ios::binary );
        if ( !input ) throw std::runtime_error( "cannot open dump input: " + path );
        std::array<char, 8> magic{};
        input.read( magic.data(), magic.size() );
        if ( magic != oneflow_1d::Weno5DumpMagic ) throw std::runtime_error( "invalid WENO5 dump magic" );
        std::uint32_t version = 0, caseCount = 0;
        oneflow_1d::ReadPod( input, version ); oneflow_1d::ReadPod( input, caseCount );
        if ( version != oneflow_1d::Weno5DumpVersion || caseCount != 2 ) throw std::runtime_error( "unsupported WENO5 dump" );
        bool pass = true;
        for ( std::uint32_t index = 0; index < caseCount; ++ index )
        {
            oneflow_1d::Weno5CaseHeader header;
            oneflow_1d::ReadPod( input, header );
            pass = RunCase( input, header ) && pass;
        }
        std::cout << "OneFLOW 1D Lax-WENO5 CPU/HIP replay: " << ( pass ? "PASS" : "FAIL" ) << "\n";
        return pass ? 0 : 1;
    }
    catch ( const std::exception & error )
    {
        std::cerr << "WENO5 HIP replay failed: " << error.what() << "\n";
        return 1;
    }
}
