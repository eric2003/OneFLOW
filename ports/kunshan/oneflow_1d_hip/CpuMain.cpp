#include "OneDUpwind.h"
#include "TestData.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    constexpr int nx = 4097;
    constexpr int uStart = -1;
    constexpr double dx = 1.0 / nx;
    const std::vector<double> velocities{ 1.0, -1.0, 0.35, -0.75 };
    const std::vector<double> u = MakeOneDState( nx, uStart );
    double maxError = 0.0;

    for ( double c : velocities )
    {
        std::vector<double> left( nx + 1 );
        std::vector<double> right( nx + 1 );
        std::vector<double> residual( nx );
        OneDCpuUpwindResidual(
            u.data(), static_cast<int>( u.size() ), uStart, nx, c, dx,
            left.data(), right.data(), residual.data() );

        const double cp = 0.5 * ( c + std::abs( c ) );
        const double cn = 0.5 * ( c - std::abs( c ) );
        for ( int i = 0; i <= nx; ++ i )
        {
            maxError = std::max(
                maxError, std::abs( left[ i ] - u[ i - 1 - uStart ] ) );
            maxError = std::max(
                maxError, std::abs( right[ i ] - u[ i - uStart ] ) );
        }
        for ( int i = 0; i < nx; ++ i )
        {
            const double expected = -(
                cp * u[ i - uStart ] + cn * u[ i + 1 - uStart ]
                - cp * u[ i - 1 - uStart ] - cn * u[ i - uStart ] ) / dx;
            maxError = std::max( maxError, std::abs( residual[ i ] - expected ) );
        }
    }

    std::cout << "OneFLOW 1D Upwind1 CPU self-test: "
              << ( maxError <= 1.0e-15 ? "PASS" : "FAIL" )
              << ", max error = " << maxError << '\n';
    return maxError <= 1.0e-15 ? 0 : 1;
}
