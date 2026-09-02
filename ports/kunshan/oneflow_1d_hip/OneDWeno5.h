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

struct Weno5HipState;

Weno5HipState * CreateWeno5HipState(
    const double * state, int nx, double gamma, double dx,
    EulerBoundary boundary );

void DestroyWeno5HipState( Weno5HipState * state ) noexcept;

void ResetWeno5HipState( Weno5HipState & state, const double * values );

void StepWeno5HipState(
    Weno5HipState & state, double dt, Weno5Trace * trace = nullptr );

void DownloadWeno5HipState( const Weno5HipState & state, double * values );

void OneDHipLaxWeno5Step(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, Weno5Trace & trace );

} // namespace oneflow_1d
