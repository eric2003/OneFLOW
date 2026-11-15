#pragma once

#include "AccelViews.h"
// The primitive-face conversion and pack are execution-space independent and
// now live in their own header, so the accelerator backend can reach them
// without including a header named after the CPU adapter.
#include "PrimitiveFaceContract.h"

BeginNameSpace( ONEFLOW )

class FluxBackend;

class EulerCpuAdapter
{
public:
    // The primitive-to-conserved conversion and equation-major pack are
    // execution-space independent. Production callers inject CPU or HIP;
    // keeping this seam here prevents the two paths from drifting.
    void CalcInvFlux(
        const PrimitiveFaceStateView & state,
        FaceFluxView & flux,
        FluxBackend & backend,
        int scheme = 0 ) const;

    void CalcInvFlux(
        const PrimitiveFaceStateView & state,
        FaceFluxView & flux,
        int scheme = 0 ) const;

    void AddFaceFlux(
        const FaceFluxView & flux,
        const FaceConnectivityView & connectivity,
        ResidualView & residual ) const;
};

EndNameSpace
