#pragma once
#include <cmath>
#include <vector>

inline std::vector<double> MakeOneDState( int nx, int uStart )
{
    const int uEnd = nx;
    std::vector<double> u( uEnd - uStart + 1 );
    for ( int i = uStart; i <= uEnd; ++ i )
    {
        const double x = static_cast<double>( i ) / static_cast<double>( nx );
        u[ i - uStart ] = 0.75 + 0.2 * std::sin( 6.283185307179586 * x )
            + 0.05 * std::cos( 18.84955592153876 * x );
    }
    return u;
}
