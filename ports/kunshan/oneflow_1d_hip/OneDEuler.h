#pragma once

#include <vector>

namespace oneflow_1d
{

constexpr int EulerComponents = 3;
constexpr int EulerRkStages = 3;

enum class EulerBoundary
{
    Periodic = 0,
    Transmissive = 1
};

struct EulerTrace
{
    int nx = 0;
    std::vector<double> faceLeft;
    std::vector<double> faceRight;
    std::vector<double> numericalFlux;
    std::vector<double> residual;
    std::vector<double> state;
};

void ResizeEulerTrace( EulerTrace & trace, int nx );

void OneDCpuEulerStep(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace );

void OneDHipEulerStep(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace );

} // namespace oneflow_1d
