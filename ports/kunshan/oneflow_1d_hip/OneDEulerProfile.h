#pragma once

#include "OneDEuler.h"

namespace oneflow_1d
{

struct EulerProfileTiming
{
    double hostMilliseconds = 0.0;
    double allocationMilliseconds = 0.0;
    double h2dMilliseconds = 0.0;
    double kernelMilliseconds = 0.0;
    double d2hMilliseconds = 0.0;
    int allocationCount = 0;
    int kernelLaunches = 0;
    bool deviceTimingAvailable = false;
    char deviceName[ 256 ] {};
    char architecture[ 64 ] {};
};

void ProfiledHipEulerStep(
    const double * state, int nx, double gamma, double dt, double dx,
    EulerBoundary boundary, EulerTrace & trace, EulerProfileTiming & timing );

} // namespace oneflow_1d
