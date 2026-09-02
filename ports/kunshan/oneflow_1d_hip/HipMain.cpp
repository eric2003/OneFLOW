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
    double maxLeftError = 0.0;
    double maxRightError = 0.0;
    double maxResidualError = 0.0;

    for ( double c : velocities )
    {
        std::vector<double> cpuLeft( nx + 1 );
        std::vector<double> cpuRight( nx + 1 );
        std::vector<double> cpuResidual( nx );
        std::vector<double> hipLeft( nx + 1 );
        std::vector<double> hipRight( nx + 1 );
        std::vector<double> hipResidual( nx );
        OneDCpuUpwindResidual(
            u.data(), static_cast<int>( u.size() ), uStart, nx, c, dx,
            cpuLeft.data(), cpuRight.data(), cpuResidual.data() );
        OneDHipUpwindResidual(
            u.data(), static_cast<int>( u.size() ), uStart, nx, c, dx,
            hipLeft.data(), hipRight.data(), hipResidual.data() );

        for ( int i = 0; i <= nx; ++ i )
        {
            maxLeftError = std::max(
                maxLeftError, std::abs( cpuLeft[ i ] - hipLeft[ i ] ) );
            maxRightError = std::max(
                maxRightError, std::abs( cpuRight[ i ] - hipRight[ i ] ) );
        }
        for ( int i = 0; i < nx; ++ i )
        {
            maxResidualError = std::max(
                maxResidualError,
                std::abs( cpuResidual[ i ] - hipResidual[ i ] ) );
        }
    }

    const double maxError = std::max(
        maxLeftError, std::max( maxRightError, maxResidualError ) );
    std::cout << "OneFLOW 1D Upwind1 HIP validation: "
              << ( maxError <= 1.0e-15 ? "PASS" : "FAIL" )
              << ", left = " << maxLeftError
              << ", right = " << maxRightError
              << ", residual = " << maxResidualError << '\n';
    return maxError <= 1.0e-15 ? 0 : 1;
}
