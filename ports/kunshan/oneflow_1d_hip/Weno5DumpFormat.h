#pragma once

#include "OneDWeno5.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace oneflow_1d
{

constexpr std::array<char, 8> Weno5DumpMagic{{'O','F','W','5','D','M','P','1'}};
constexpr std::uint32_t Weno5DumpVersion = 1;

struct Weno5CaseHeader
{
    std::array<char, 32> name{};
    std::int32_t nx = 0;
    std::int32_t steps = 0;
    std::int32_t boundary = 0;
    std::int32_t reserved = 0;
    double gamma = 0.0;
    double dt = 0.0;
    double dx = 0.0;
};

template<class T>
void WritePod( std::ostream & output, const T & value )
{
    output.write( reinterpret_cast<const char *>( &value ), sizeof( T ) );
    if ( ! output ) throw std::runtime_error( "failed to write WENO5 dump" );
}

template<class T>
void ReadPod( std::istream & input, T & value )
{
    input.read( reinterpret_cast<char *>( &value ), sizeof( T ) );
    if ( ! input ) throw std::runtime_error( "failed to read WENO5 dump" );
}

inline void WriteValues( std::ostream & output, const std::vector<double> & values )
{
    output.write( reinterpret_cast<const char *>( values.data() ),
                  static_cast<std::streamsize>( values.size() * sizeof( double ) ) );
    if ( ! output ) throw std::runtime_error( "failed to write WENO5 values" );
}

inline void ReadValues( std::istream & input, std::vector<double> & values )
{
    input.read( reinterpret_cast<char *>( values.data() ),
                static_cast<std::streamsize>( values.size() * sizeof( double ) ) );
    if ( ! input ) throw std::runtime_error( "failed to read WENO5 values" );
}

inline void WriteTrace( std::ostream & output, const Weno5Trace & trace )
{
    WriteValues( output, trace.splitPositive );
    WriteValues( output, trace.splitNegative );
    WriteValues( output, trace.reconstructedPositive );
    WriteValues( output, trace.reconstructedNegative );
    WriteValues( output, trace.numericalFlux );
    WriteValues( output, trace.residual );
    WriteValues( output, trace.state );
}

inline void ReadTrace( std::istream & input, int nx, Weno5Trace & trace )
{
    ResizeWeno5Trace( trace, nx );
    ReadValues( input, trace.splitPositive );
    ReadValues( input, trace.splitNegative );
    ReadValues( input, trace.reconstructedPositive );
    ReadValues( input, trace.reconstructedNegative );
    ReadValues( input, trace.numericalFlux );
    ReadValues( input, trace.residual );
    ReadValues( input, trace.state );
}

} // namespace oneflow_1d
