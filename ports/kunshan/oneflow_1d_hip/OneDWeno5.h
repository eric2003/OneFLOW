#pragma once

#include "OneDEuler.h"

#include <vector>

namespace oneflow_1d
{

constexpr int WenoGhostCells = 3;

struct Weno5Trace
{
    int nx = 0;
    std::vector<double> splitPositive;
    std::vector<double> splitNegative;
    std::vector<double> reconstructedPositive;
    std::vector<double> reconstructedNegative;
    std::vector<double> numericalFlux;
    std::vector<double> residual;
    std::vector<double> state;
};

void ResizeWeno5Trace( Weno5Trace & trace, int nx );

void OneDCpuLaxWeno5Step(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, Weno5Trace & trace );

void OneDHipLaxWeno5Step(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, Weno5Trace & trace );

} // namespace oneflow_1d
