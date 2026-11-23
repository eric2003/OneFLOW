/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// Single-faceside Barth limiter contract between the device kernel and the
// host policy functions that the CPU oracle uses.
//
// The limiter is the first capability that keeps the flagship case off the
// device, so its formula must be pinned: the phi produced for one face side
// has to match BarthFunction(ReconstructFaceIncrement(...), ..., FaceSideDirectionDot(...))
// including the sign convention and the three branches (dot < 0, dq above the
// range, dq below it).
//
// Scope: one face side.  The mesh-level limiter additionally needs the cell
// neighbourhood range and the minimum over faces; those are covered by the
// cross-backend gate once the limiter is wired into reconstruction.

#include "AccelRuntime.h"
#include "HipFluxBackend.h"
#include "ULimiterPolicy.h"

#include <gtest/gtest.h>

#include <cmath>

using namespace ONEFLOW;

namespace
{

struct SideSample
{
    double dqdx, dqdy, dqdz;
    double dx, dy, dz;
    double nx, ny, nz;
    double minv, maxv;
    double orientationSign;
    const char * name;
};

constexpr double kTolerance = 1.0e-14;

class DeviceLimiterContract : public ::testing::Test
{
protected:
    static void SetUpTestSuite()
    {
        InitializeAccelRuntime( 0, 1 );
    }

    static void TearDownTestSuite()
    {
        AccelRuntime::Instance().Finalize();
    }

    static void Check( const SideSample & sample )
    {
        const double hostIncrement = ReconstructFaceIncrement(
            sample.dqdx, sample.dqdy, sample.dqdz,
            sample.dx, sample.dy, sample.dz );
        const double hostDot = FaceSideDirectionDot(
            sample.nx, sample.ny, sample.nz,
            sample.dx, sample.dy, sample.dz,
            sample.orientationSign );
        const double hostPhi = BarthFunction(
            hostIncrement, sample.minv, sample.maxv, hostDot );

        const double devicePhi = HipFluxBackend::Shared().EvaluateDeviceBarth(
            sample.dqdx, sample.dqdy, sample.dqdz,
            sample.dx, sample.dy, sample.dz,
            sample.nx, sample.ny, sample.nz,
            sample.minv, sample.maxv, sample.orientationSign );

        EXPECT_NEAR( devicePhi, hostPhi, kTolerance ) << sample.name;
    }
};

}

TEST_F( DeviceLimiterContract, InsideRangeKeepsFullFactor )
{
    const SideSample sample{
        0.1, -0.05, 0.02, 0.5, 0.25, -0.1, 1.0, 0.0, 0.0,
        -0.05, 0.05, 1.0, "inside-range" };
    Check( sample );
}

TEST_F( DeviceLimiterContract, AboveRangeTightens )
{
    const SideSample sample{
        1.0, 0.4, 0.2, 0.5, 0.25, 0.1, 1.0, 0.0, 0.0,
        -0.05, 0.10, 1.0, "above-range" };
    Check( sample );
}

TEST_F( DeviceLimiterContract, BelowRangeTightens )
{
    const SideSample sample{
        -0.8, -0.3, 0.1, 0.5, 0.25, 0.1, 1.0, 0.0, 0.0,
        -0.10, 0.05, 1.0, "below-range" };
    Check( sample );
}

TEST_F( DeviceLimiterContract, OppositeOrientationCollapses )
{
    // dot < 0 returns zero regardless of the increment.
    const SideSample sample{
        0.1, 0.0, 0.0, -0.5, 0.0, 0.0, 1.0, 0.0, 0.0,
        -0.05, 0.05, 1.0, "opposite-orientation" };
    Check( sample );
}

TEST_F( DeviceLimiterContract, SkewedFaceUsesNormalizedDirection )
{
    const SideSample sample{
        0.3, -0.2, 0.15, 0.2, -0.4, 0.3, 0.6, 0.5, -0.3,
        -0.02, 0.04, -1.0, "skewed-face" };
    Check( sample );
}
